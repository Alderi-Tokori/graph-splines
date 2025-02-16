#[derive(Debug, Copy, Clone)]
pub struct Point {
    x: f64,
    y: f64,
}

trait PolynomialFunction {
    fn eval(&self, x: f64) -> Option<f64>;
}

#[derive(Debug)]
struct Polynomial {
    coefficients: Vec<f64>,
}

impl PolynomialFunction for Polynomial {
    fn eval(&self, x: f64) -> Option<f64> {
        if self.coefficients.len() == 0 {
            return None;
        }

        let mut sum = 0.0;
        let mut x_degree = 1.0;

        self.coefficients.iter().for_each(|c| {
            sum += c * x_degree;
            x_degree *= x;
        });

        Some(sum)
    }
}

#[derive(Debug)]
struct GraphSplineInterval {
    start: Point,
    end: Point,
    polynomial: Polynomial,
}

impl PolynomialFunction for GraphSplineInterval {
    fn eval(&self, x: f64) -> Option<f64> {
        self.polynomial.eval(x - self.start.x)
    }
}

#[derive(Debug)]
struct GraphSpline {
    intervals: Vec<GraphSplineInterval>
}

impl PolynomialFunction for GraphSpline {
    fn eval(&self, x: f64) -> Option<f64> {
        if self.intervals.len() == 0 {
            return None
        }

        let first_interval = self.intervals.first().unwrap();
        let last_interval = self.intervals.last().unwrap();

        if x < first_interval.end.x {
            return first_interval.eval(x);
        } else if x >= last_interval.start.x {
            return last_interval.eval(x);
        }

        self.intervals[1..self.intervals.len() - 1]
            .into_iter()
            .find(|i| {
                i.start.x <= x && i.end.x > x
            })
            .expect("We've tested the two edge cases already so we should always have an interval here, as all intervals should be next to each other")
            .eval(x)
    }
}

fn get_graph_spline_intervals(points: &Vec<Point>) -> Vec<GraphSplineInterval> {
    let mut interval_degrees_list= Vec::new();

    points
        .windows(3)
        .for_each(|points| {
            let mut degree = 3;

            if (points[0].y < points[1].y && points[2].y <= points[1].y)
                || (points[0].y > points[1].y && points[2].y >= points[1].y) {
                // To accomodate for the additionnal constraint at local optimum, the splines leading to a local
                // optimum will need to be quartic instead of only cubic
                degree += 1;
            }

            interval_degrees_list.push(GraphSplineInterval {
                start: points[0],
                end: points[1],
                polynomial: Polynomial {
                    coefficients: vec![0.0; degree + 1]
                }
            })
        })
    ;

    let nb_points = points.len();
    if nb_points >= 2 {
        interval_degrees_list.push(GraphSplineInterval {
            start: points[points.len() - 2],
            end: points[points.len() - 1],
            polynomial: Polynomial {
                coefficients: vec![0.0; 4]
            }
        })
    }

    interval_degrees_list
}

fn get_equation_factors(equation: &mut Vec<f64>, coefficient_idx: &usize, nb_coefficients: &usize, x_value: &f64, derivative: &usize) {
    let mut cur_x_value = 1.0;
    for i in *derivative..*nb_coefficients {
        let mut derivative_coeff = 1.0;
        for j in 0..*derivative {
            derivative_coeff = derivative_coeff * ((i - j) as f64);
        }

        equation[coefficient_idx + i] = 1.0 * derivative_coeff * cur_x_value;

        cur_x_value *= x_value;
    }
}

fn get_graph_spline_equation_matrix(intervals: &Vec<GraphSplineInterval>) -> Vec<Vec<f64>> {
    let nb_unknowns = intervals
        .iter()
        .map(|i| i.polynomial.coefficients.len())
        .sum()
    ;

    let mut equation_matrix = vec![vec![0.0; nb_unknowns + 1]; nb_unknowns];
    let mut equation_idx = 0;
    let mut coefficient_idx = 0;

    // Initial boundary condition
    // Natural spline
    equation_matrix[equation_idx][coefficient_idx + 2] = 2.0;
    equation_idx += 1;
    
    // Starts on first point
    let first_interval = intervals.first().expect("We should have at least one interval if we called this function");
    equation_matrix[equation_idx][coefficient_idx] = 1.0;
    equation_matrix[equation_idx][nb_unknowns] = first_interval.start.y;
    equation_idx += 1;

    // inner knots conditions
    intervals
        .windows(2)
        .for_each(|intervals| {
            let max_x_value_left = intervals[0].end.x - intervals[0].start.x;
            let coefficient_idx_right = coefficient_idx + intervals[0].polynomial.coefficients.len();

            // C0 continuity left interval
            get_equation_factors(
                &mut equation_matrix[equation_idx],
                &coefficient_idx,
                &intervals[0].polynomial.coefficients.len(),
                &max_x_value_left,
                &0
            );
            equation_matrix[equation_idx][nb_unknowns] = intervals[0].end.y;
            equation_idx += 1;
            
            // C0 continuity right interval
            get_equation_factors(
                &mut equation_matrix[equation_idx],
                &coefficient_idx_right,
                &intervals[1].polynomial.coefficients.len(),
                &0.0,
                &0
            );
            equation_matrix[equation_idx][nb_unknowns] = intervals[1].start.y;
            equation_idx += 1;
            
            // C1 continuity 
            let derivative_number = 1;
            get_equation_factors(
                &mut equation_matrix[equation_idx],
                &coefficient_idx,
                &intervals[0].polynomial.coefficients.len(),
                &max_x_value_left,
                &derivative_number
            );
            get_equation_factors(
                &mut equation_matrix[equation_idx],
                &coefficient_idx_right,
                &intervals[1].polynomial.coefficients.len(),
                &0.0,
                &derivative_number
            );
            equation_idx += 1;
            
            // C2 continuity
            let derivative_number = 2;
            get_equation_factors(
                &mut equation_matrix[equation_idx],
                &coefficient_idx,
                &intervals[0].polynomial.coefficients.len(),
                &max_x_value_left,
                &derivative_number
            );
            get_equation_factors(
                &mut equation_matrix[equation_idx],
                &coefficient_idx_right,
                &intervals[1].polynomial.coefficients.len(),
                &0.0,
                &derivative_number
            );
            equation_idx += 1;
            
            // Additional constraint for local optimum: first derivative must be zero
            if intervals[0].polynomial.coefficients.len() == 5 {
                get_equation_factors(
                    &mut equation_matrix[equation_idx],
                    &coefficient_idx,
                    &intervals[0].polynomial.coefficients.len(),
                    &max_x_value_left,
                    &1
                );
                equation_idx += 1;
            }

            coefficient_idx += intervals[0].polynomial.coefficients.len();
        })
    ;

    // Final boundary conditions
    // Natural spline
    let last_interval = intervals.last().expect("We should have at least one interval if we called this function");
    let max_x_value = last_interval.end.x - last_interval.start.x;
    get_equation_factors(
        &mut equation_matrix[equation_idx],
        &coefficient_idx,
        &last_interval.polynomial.coefficients.len(),
        &max_x_value,
        &2
    );
    equation_idx += 1;

    // Ends on last point
    get_equation_factors(
        &mut equation_matrix[equation_idx],
        &coefficient_idx,
        &last_interval.polynomial.coefficients.len(),
        &max_x_value,
        &0
    );
    equation_matrix[equation_idx][nb_unknowns] = last_interval.end.y;

    equation_matrix
}

// pub fn get_graph_spline_interpolation_function(points: Vec<Point>) -> Option<GraphSpline> {
//     if points.len() < 2 {
//         return None;
//     }
//     
//     // We're going to construct the equation system from the points
//     let intervals = get_graph_spline_intervals(&points);
//     let equation_matrix = get_graph_spline_equation_matrix(&intervals);
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_correctly_evaluates_polynomial() {
        let p = Polynomial {
            coefficients: vec![5.0, 1.0, 1.0]
        };

        assert_eq!(p.eval(0.0), Some(5.0));
        assert_eq!(p.eval(10.0), Some(115.0));
    }

    #[test]
    fn it_correctly_returns_none_for_undefined_polynomial() {
        let p = Polynomial {
            coefficients: vec![]
        };

        assert_eq!(p.eval(0.0), None);
    }

    #[test]
    fn it_correctly_evaluates_piecewise_polynomial() {
        let p = GraphSpline {
            intervals: vec![
                GraphSplineInterval {
                    start: Point {x: 0.0, y: 0.0},
                    end: Point {x: 5.0, y: 0.0},
                    polynomial: Polynomial {
                        coefficients: vec![5.0, 2.0]
                    }
                },
                GraphSplineInterval {
                    start: Point {x: 5.0, y: 0.0},
                    end: Point {x: 10.0, y: 0.0},
                    polynomial: Polynomial {
                        coefficients: vec![0.0, 1.0]
                    }
                },
                GraphSplineInterval {
                    start: Point {x: 10.0, y: 0.0},
                    end: Point {x: 15.0, y: 0.0},
                    polynomial: Polynomial {
                        coefficients: vec![0.0, 1.0, 1.0]
                    }
                }
            ]
        };

        assert_eq!(p.eval(0.0), Some(5.0));
        assert_eq!(p.eval(4.99), Some(14.98));
        assert_eq!(p.eval(5.0), Some(0.0));
        assert_eq!(p.eval(9.99), Some(4.99));
        assert_eq!(p.eval(10.0), Some(0.0));
        assert_eq!(p.eval(12.0), Some(6.0));
        assert_eq!(p.eval(-1.0), Some(3.0));
        assert_eq!(p.eval(20.0), Some(110.0));
    }

    #[test]
    fn it_correctly_returns_none_if_no_interval_defined() {
        let p = GraphSpline {
            intervals: vec![]
        };

        assert_eq!(p.eval(1.0), None);
    }

    #[test]
    fn it_correctly_calculates_interval_degrees_list() {
        let points = vec![
            Point {x: 0.0, y: 10.0},
            Point {x: 1.0, y: 20.0},
            Point {x: 2.0, y: 30.0},
            Point {x: 3.0, y: 20.0},
            Point {x: 4.0, y: 10.0},
            Point {x: 5.0, y: 20.0},
            Point {x: 6.0, y: 20.0},
            Point {x: 7.0, y: 10.0},
            Point {x: 8.0, y: 0.0},
            Point {x: 9.0, y: 0.0},
            Point {x: 10.0, y: 10.0},
        ];

        let intervals = get_graph_spline_intervals(&points);

        assert_eq!(intervals[0].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[1].polynomial.coefficients.len(), 5);
        assert_eq!(intervals[2].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[3].polynomial.coefficients.len(), 5);
        assert_eq!(intervals[4].polynomial.coefficients.len(), 5);
        assert_eq!(intervals[5].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[6].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[7].polynomial.coefficients.len(), 5);
        assert_eq!(intervals[8].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[9].polynomial.coefficients.len(), 4);
    }

    #[test]
    fn it_correctly_generates_equation_matrix() {
        let points = vec![
            Point {x: 0.0, y: 1.0},
            Point {x: 1.0, y: 3.0},
            Point {x: 2.0, y: 2.0}
        ];
        
        let intervals = get_graph_spline_intervals(&points);
        
        let equation_matrix = get_graph_spline_equation_matrix(&intervals);
        assert_eq!(equation_matrix, vec![
            vec![0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            vec![1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 3.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 3.0],
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 2.0, 6.0, 12.0, 0.0, 0.0, 2.0, 0.0, 0.0],
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 6.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0],
        ])
    }
}
