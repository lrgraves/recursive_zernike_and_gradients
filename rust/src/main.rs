use nalgebra as na;
use std::collections::HashMap;

fn main() {
    // Define an example of input coordinates
    let coordinates = vec![
        na::Point2::new(0.0, 1.0),
        na::Point2::new(0.0, 0.5),
        na::Point2::new(0.0, 0.0),
        na::Point2::new(0.0, -0.5),
        na::Point2::new(0.0, -1.0),
    ];

    // Define an example maximum radial order, here 4, which corresponds to
    // a maximum ansi index for the zernike calculated of 14
    let max_order = 4;

    // Then, calculate our Zernike function, and gradients of it.
    let (mut zernike_values, mut gradient_x_values, mut gradient_y_values) =
        gradient(&coordinates, max_order);

    // The output is not ordered, I am sure there is some clever rust way to avoid this initially, but I don't know it, thus we sort
    zernike_values.sort_by(|a, b| Ord::cmp(&a.0.j, &b.0.j));
    gradient_x_values.sort_by(|a, b| Ord::cmp(&a.0.j, &b.0.j));
    gradient_y_values.sort_by(|a, b| Ord::cmp(&a.0.j, &b.0.j));

    // And lastly, print out our results (if you want, there is no real purpose for printing it out however besides cluttering your terminal)
    print!("The zernike val, dx, dy is:");
    for index in 0..zernike_values.len() {
        for pos in 0..coordinates.len() {
            print!(
                "z: {}, dz/dx: {}, dz/dy: {} at (x: {}, y: {}) for ansi index: {}, ",
                zernike_values[index].1[pos],
                gradient_x_values[index].1[pos],
                gradient_y_values[index].1[pos],
                coordinates[pos].x,
                coordinates[pos].y,
                zernike_values[index].0.j,
            );
        }
    }
}

#[derive(Debug, Clone)]
pub struct Ansi {
    j: u16,
}

pub fn gradient(
    coordinates: &[na::Point2<f64>],
    mut max_order: u16,
) -> (
    Vec<(Ansi, Vec<f64>)>,
    Vec<(Ansi, Vec<f64>)>,
    Vec<(Ansi, Vec<f64>)>,
) {
    //   Calculates recursively the Zernike and Gradient polynomialt sets
    //  The zernike polynomial gradients are ouput as a vector of the ansi order
    // and the evaluated gradient value at the provided input coordinates
    // This algorithm is taken from the following publication:
    // Torben B. Andersen, "Efficient and robust recurrence relations for the Zernike circle polynomials and their derivatives in Cartesian coordinates,"
    // Opt. Express 26, 18878-18896 (2018)
    // Author : Logan Rodriguez Graves. Date: 1/27/022

    //
    // Inputs:
    // coordinates: (x,y) point coordinates of where you would like to evaluate your Zernike, and Gradients functions
    // max_order: maximum radial order to recursively generate the Zernike, and Zernike Gradient Functions to.
    // Outputs a tuple of vectors of tuples: (Vec<(Ansi Order, Vec<zernike values>)>, Vec<(Ansi Order, Vec<x gradient values>)>, Vec<(Ansi Order, Vec<y gradient values>)>)
    // where each function has been evaluated at the input coordinate points for each ansi order

    // Determine the max radial order, if it is less than 2 then set it to at least 2.
    if max_order < 2 {
        max_order = 2;
    }

    // Hashmaps are your friend here! They will allow us to directly reference into
    // the sets we need, as opposed to an unnecessarily extra complex step of translating nm
    // indices to a single order.
    let mut zernike_map: HashMap<(u16, u16), (Ansi, Vec<f64>)> = HashMap::new();
    let mut gradient_x_map: HashMap<(u16, u16), (Ansi, Vec<f64>)> = HashMap::new();
    let mut gradient_y_map: HashMap<(u16, u16), (Ansi, Vec<f64>)> = HashMap::new();

    // Define the seeding Ansi order 1 and 2 zernike polynomials
    let mut seed_zernike_1 = vec![1.0; coordinates.len()];
    let mut seed_zernike_2 = vec![1.0; coordinates.len()];

    for pos in 0..coordinates.len() {
        seed_zernike_1[pos] *= coordinates[pos].y;
        seed_zernike_2[pos] *= coordinates[pos].x;
    }

    // Define a variable for the ansi order going forward
    let mut ansi_order = 2;

    // Add in the seed values to our hashmaps
    zernike_map.insert((0, 0), (Ansi { j: 0 }, vec![1.0; coordinates.len()]));
    zernike_map.insert((1, 0), (Ansi { j: 1 }, seed_zernike_1));
    zernike_map.insert((1, 1), (Ansi { j: 2 }, seed_zernike_2));

    gradient_x_map.insert((0, 0), (Ansi { j: 0 }, vec![0.0; coordinates.len()]));
    gradient_x_map.insert((1, 0), (Ansi { j: 1 }, vec![0.0; coordinates.len()]));
    gradient_x_map.insert((1, 1), (Ansi { j: 2 }, vec![1.0; coordinates.len()]));

    gradient_y_map.insert((0, 0), (Ansi { j: 0 }, vec![0.0; coordinates.len()]));
    gradient_y_map.insert((1, 0), (Ansi { j: 1 }, vec![1.0; coordinates.len()]));
    gradient_y_map.insert((1, 1), (Ansi { j: 2 }, vec![0.0; coordinates.len()]));

    // Outer loop for radial order index
    for n in 2..max_order + 1 {
        // Inner loop for azimuthal index
        for m in 0..n + 1 {
            // Step the ansi order for our next iterations
            ansi_order += 1;
            // Define the temporary vectors to hold the new calculated zernike and gradient values
            // They get consumed in each m loop into a hashmap, hence instantiating them here
            let mut new_zern_values = Vec::with_capacity(coordinates.len());
            let mut new_gradient_x = Vec::with_capacity(coordinates.len());
            let mut new_gradient_y = Vec::with_capacity(coordinates.len());

            // Handle the various cases of exception
            if m == 0 {
                // Collect Needed Prior Zernike and Gradient Polynomials
                // Index (n - 1, 0)
                let (_, zern_n1_0) = zernike_map.get(&(n - 1, 0)).unwrap();
                // Index (n - 1, n - 1)
                let (_, zern_n1_n1) = zernike_map.get(&(n - 1, n - 1)).unwrap();

                // Iterating through every input coordinate point, define the new
                // Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in 0..coordinates.len() {
                    let x: f64 = coordinates[pos].x;
                    let y: f64 = coordinates[pos].y;

                    new_zern_values.push(x * zern_n1_0[pos] + y * zern_n1_n1[pos]);
                    new_gradient_x.push(n as f64 * zern_n1_0[pos]);
                    new_gradient_y.push(n as f64 * zern_n1_n1[pos]);
                }
            } else if m == n {
                // Collect Needed Prior Zernike and Gradient Polynomials
                // Index (n - 1, 0)
                let (_, zern_n1_0) = zernike_map.get(&(n - 1, 0)).unwrap();
                // Index (n - 1, n - 1)
                let (_, zern_n1_n1) = zernike_map.get(&(n - 1, n - 1)).unwrap();

                // Iterating through every input coordinate point, define the new
                // Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in 0..coordinates.len() {
                    let x: f64 = coordinates[pos].x;
                    let y: f64 = coordinates[pos].y;

                    new_zern_values.push(x * zern_n1_n1[pos] - y * zern_n1_0[pos]);
                    new_gradient_x.push(n as f64 * zern_n1_n1[pos]);
                    new_gradient_y.push(-1.0 * n as f64 * zern_n1_0[pos]);
                }
            } else if n % 2 != 0 && m == (n - 1) / 2 {
                // Collect Needed Prior Zernike and Gradient Polynomials
                // Index (n - 1, n - 1 - m)
                let (_, zern_n1_n1m) = zernike_map.get(&(n - 1, n - 1 - m)).unwrap();

                // Index (n - 1, m - 1)
                let (_, zern_n1_m1) = zernike_map.get(&(n - 1, m - 1)).unwrap();

                // Index (n - 1, n - m)
                let (_, zern_n1_nm) = zernike_map.get(&(n - 1, n - m)).unwrap();

                // Index (n - 2, m - 1)
                let (_, zern_n2_m1) = zernike_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_x_n2_m1) = gradient_x_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_y_n2_m1) = gradient_y_map.get(&(n - 2, m - 1)).unwrap();

                // Iterating through every input coordinate point, define the new
                // Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in 0..coordinates.len() {
                    let x: f64 = coordinates[pos].x;
                    let y: f64 = coordinates[pos].y;

                    // Calculate new values
                    new_zern_values.push(
                        y * zern_n1_n1m[pos] + x * zern_n1_m1[pos]
                            - y * zern_n1_nm[pos]
                            - zern_n2_m1[pos],
                    );

                    new_gradient_x.push(n as f64 * zern_n1_m1[pos] + gradient_x_n2_m1[pos]);
                    new_gradient_y.push(
                        n as f64 * zern_n1_n1m[pos] - n as f64 * zern_n1_nm[pos]
                            + gradient_y_n2_m1[pos],
                    );
                }
            } else if n % 2 != 0 && m == (n - 1) / 2 + 1 {
                // Collect Needed Prior Zernike and Gradient Polynomials

                // Index (n - 1, n - 1 - m)
                let (_, zern_n1_n1m) = zernike_map.get(&(n - 1, n - 1 - m)).unwrap();

                // Index (n - 1, m)
                let (_, zern_n1_m) = zernike_map.get(&(n - 1, m)).unwrap();

                // Index (n - 1, m - 1)
                let (_, zern_n1_m1) = zernike_map.get(&(n - 1, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, zern_n2_m1) = zernike_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_x_n2_m1) = gradient_x_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_y_n2_m1) = gradient_y_map.get(&(n - 2, m - 1)).unwrap();

                // Iterating through every input coordinate point, define the new
                // Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in 0..coordinates.len() {
                    let x: f64 = coordinates[pos].x;
                    let y: f64 = coordinates[pos].y;

                    // Calculate new values
                    new_zern_values.push(
                        x * zern_n1_m[pos] + y * zern_n1_n1m[pos] + x * zern_n1_m1[pos]
                            - zern_n2_m1[pos],
                    );

                    new_gradient_x.push(
                        n as f64 * zern_n1_m[pos]
                            + n as f64 * zern_n1_m1[pos]
                            + gradient_x_n2_m1[pos],
                    );
                    new_gradient_y.push(n as f64 * zern_n1_n1m[pos] + gradient_y_n2_m1[pos]);
                }
            } else if n % 2 == 0 && m == n / 2 {
                // Collect Needed Prior Zernike and Gradient Polynomials

                // Index (n - 1, n - 1 - m)
                let (_, zern_n1_n1m) = zernike_map.get(&(n - 1, n - 1 - m)).unwrap();

                // Index (n - 1, m)
                let (_, zern_n1_m) = zernike_map.get(&(n - 1, m)).unwrap();

                // Index (n - 1, m - 1)
                let (_, zern_n1_m1) = zernike_map.get(&(n - 1, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, zern_n2_m1) = zernike_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_x_n2_m1) = gradient_x_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_y_n2_m1) = gradient_y_map.get(&(n - 2, m - 1)).unwrap();

                // Iterating through every input coordinate point, define the new
                // Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in 0..coordinates.len() {
                    let x: f64 = coordinates[pos].x;
                    let y: f64 = coordinates[pos].y;

                    // Calculate new values
                    new_zern_values.push(
                        2.0 * x * zern_n1_m[pos] + 2.0 * y * zern_n1_m1[pos] - zern_n2_m1[pos],
                    );

                    new_gradient_x.push(2.0 * n as f64 * zern_n1_m[pos] + gradient_x_n2_m1[pos]);
                    new_gradient_y.push(2.0 * n as f64 * zern_n1_n1m[pos] + gradient_y_n2_m1[pos]);
                }
            } else {
                // Collect Needed Prior Zernike and Gradient Polynomials

                // Index (n - 1, n - m)
                let (_, zern_n1_nm) = zernike_map.get(&(n - 1, n - m)).unwrap();

                // Index (n - 1, n - 1 - m)
                let (_, zern_n1_n1m) = zernike_map.get(&(n - 1, n - 1 - m)).unwrap();

                // Index (n - 1, m)
                let (_, zern_n1_m) = zernike_map.get(&(n - 1, m)).unwrap();

                // Index (n - 1, m - 1)
                let (_, zern_n1_m1) = zernike_map.get(&(n - 1, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, zern_n2_m1) = zernike_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_x_n2_m1) = gradient_x_map.get(&(n - 2, m - 1)).unwrap();

                // Index (n - 2, m - 1)
                let (_, gradient_y_n2_m1) = gradient_y_map.get(&(n - 2, m - 1)).unwrap();

                // Iterating through every input coordinate point, define the new
                // Zernike, and Zernike Gradient values evaluated at the coordinates input.
                for pos in 0..coordinates.len() {
                    let x: f64 = coordinates[pos].x;
                    let y: f64 = coordinates[pos].y;

                    // Calculate new values
                    new_zern_values.push(
                        x * zern_n1_m[pos] + y * zern_n1_n1m[pos] + x * zern_n1_m1[pos]
                            - y * zern_n1_nm[pos]
                            - zern_n2_m1[pos],
                    );

                    new_gradient_x.push(
                        n as f64 * zern_n1_m[pos]
                            + n as f64 * zern_n1_m1[pos]
                            + gradient_x_n2_m1[pos],
                    );
                    new_gradient_y.push(
                        n as f64 * zern_n1_n1m[pos] - n as f64 * zern_n1_nm[pos]
                            + gradient_y_n2_m1[pos],
                    );
                }
            }
            // Add the new values to our nifty hashmap

            zernike_map.insert((n, m), (Ansi { j: ansi_order }, new_zern_values));
            gradient_x_map.insert((n, m), (Ansi { j: ansi_order }, new_gradient_x));
            gradient_y_map.insert((n, m), (Ansi { j: ansi_order }, new_gradient_y));
        } // End of inner azimuthal loop
    } // End of outer radial order loop

    // There is no real reason to convert out of a hashmap here, other than we use a vec outside of this.
    // However, I am pretty fond of the hashmap and may change it to output that instead
    let output_zern = zernike_map
        .values()
        .cloned()
        .collect::<Vec<(Ansi, Vec<f64>)>>();

    let output_gradient_x = gradient_x_map
        .values()
        .cloned()
        .collect::<Vec<(Ansi, Vec<f64>)>>();

    let output_gradient_y = gradient_y_map
        .values()
        .cloned()
        .collect::<Vec<(Ansi, Vec<f64>)>>();

    (output_zern, output_gradient_x, output_gradient_y)
}