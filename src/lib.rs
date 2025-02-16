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

// pub fn get_graph_spline_interpolation_function(points: Vec<Point>) -> GraphSpline {
//     // We're going to construct the equation system from the points
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
}
