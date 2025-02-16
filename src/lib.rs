#[derive(Debug)]
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
    start: f64,
    end: f64,
    polynomial: Polynomial,
}

impl PolynomialFunction for GraphSplineInterval {
    fn eval(&self, x: f64) -> Option<f64> {
        self.polynomial.eval(x - self.start)
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

        if x < first_interval.end {
            return first_interval.eval(x);
        } else if x >= last_interval.start {
            return last_interval.eval(x);
        }

        self.intervals[1..self.intervals.len() - 1]
            .into_iter()
            .find(|i| {
                i.start <= x && i.end > x
            })
            .expect("We've tested the two edge cases already so we should always have an interval here, as all intervals should be next to each other")
            .eval(x)
    }
}

fn get_graph_spline_intervals(points: &Vec<Point>) -> Vec<GraphSplineInterval> {
    let mut interval_degrees_list= Vec::new();

    points
        .windows(3)
        .enumerate()
        .for_each(|(idx, points)| {
            let mut degree = 3;

            if (points[0].y < points[1].y && points[2].y <= points[1].y)
                || (points[0].y > points[1].y && points[2].y >= points[1].y) {
                // To accomodate for the additionnal constraint at local optimum, the splines leading to a local
                // optimum will need to be quartic instead of only cubic
                degree += 1;
            }

            interval_degrees_list.push(GraphSplineInterval {
                start: points[0].x,
                end: points[1].x,
                polynomial: Polynomial {
                    coefficients: vec![0.0; degree]
                }
            })
        })
    ;

    let nb_points = points.len();
    if nb_points >= 2 {
        interval_degrees_list.push(GraphSplineInterval {
            start: points[points.len() - 2].x,
            end: points[points.len() - 1].x,
            polynomial: Polynomial {
                coefficients: vec![0.0; 3]
            }
        })
    }

    interval_degrees_list
}

// fn get_graph_spline_equation_matrix(points: &Vec<Point>, interval_degrees_list: &Vec<u8>) -> Vec<Vec<u64>> {
//     let nb_unknowns:u8 = interval_degrees_list
//         .iter()
//         .map(|d| {d + 1})
//         .sum()
//     ;
//
//     let mut equation_matrix = vec![vec![0; usize::from(nb_unknowns) + 1]; usize::from(nb_unknowns)];
//     let mut equation_idx = 0;
//     let mut coefficient_idx = 0;
//
//
// }

// pub fn get_graph_spline_interpolation_function(points: Vec<Point>) -> GraphSpline {
//     // We're going to construct the equation system from the points
//     let intervals = get_graph_spline_intervals(&points);
//
//     let nb_unknowns = intervals
//         .iter()
//         .map(|i| {i.polynomial.coefficients.len() + 1})
//         .sum()
//     ;
//
//     let mut equation_matrix = vec![vec![0; nb_unknowns + 1]; nb_unknowns];
//     let mut equation_idx = 0;
//     let mut coefficient_idx = 0;
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
                    start: 0.0,
                    end: 5.0,
                    polynomial: Polynomial {
                        coefficients: vec![5.0, 2.0]
                    }
                },
                GraphSplineInterval {
                    start: 5.0,
                    end: 10.0,
                    polynomial: Polynomial {
                        coefficients: vec![0.0, 1.0]
                    }
                },
                GraphSplineInterval {
                    start: 10.0,
                    end: 15.0,
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

        assert_eq!(intervals[0].polynomial.coefficients.len(), 3);
        assert_eq!(intervals[1].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[2].polynomial.coefficients.len(), 3);
        assert_eq!(intervals[3].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[4].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[5].polynomial.coefficients.len(), 3);
        assert_eq!(intervals[6].polynomial.coefficients.len(), 3);
        assert_eq!(intervals[7].polynomial.coefficients.len(), 4);
        assert_eq!(intervals[8].polynomial.coefficients.len(), 3);
        assert_eq!(intervals[9].polynomial.coefficients.len(), 3);
    }
}
