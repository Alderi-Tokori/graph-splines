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

fn calculate_interval_degrees_list(points: &Vec<Point>) -> Vec<u8> {
    // To accomodate for the additionnal constraint at local optimum, the splines leading to a local
    // optimum will need to be quartic instead of only cubic
    let mut interval_degrees_list = vec![3; points.len() - 1];

    points
        .windows(3)
        .enumerate()
        .filter(|(idx, points)|
            (points[0].y < points[1].y && points[2].y <= points[1].y)
                || (points[0].y > points[1].y && points[2].y >= points[1].y)
        )
        .for_each(|(idx, points)|
            interval_degrees_list[idx] = 4
        )
    ;

    interval_degrees_list
}

// pub fn get_graph_spline_interpolation_function(points: Vec<Point>) -> GraphSpline {
//     // We're going to construct the equation system from the points
//     let mut interval_degrees_list = calculate_interval_degrees_list(&points);
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

        let interval_degrees = calculate_interval_degrees_list(&points);

        assert_eq!(interval_degrees[0], 3);
        assert_eq!(interval_degrees[1], 4);
        assert_eq!(interval_degrees[2], 3);
        assert_eq!(interval_degrees[3], 4);
        assert_eq!(interval_degrees[4], 4);
        assert_eq!(interval_degrees[5], 3);
        assert_eq!(interval_degrees[6], 3);
        assert_eq!(interval_degrees[7], 4);
        assert_eq!(interval_degrees[8], 3);
        assert_eq!(interval_degrees[9], 3);
    }
}
