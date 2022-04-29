# recursive_zernike_gradient
Recursive Zernike Gradient Generation
Efficient and robust recurrence relations for the Zernike circle polynomials and their derivatives in Cartesian coordinates

Background
This package provides a python and rust code for generating zernike circle polynomials and their derivatives using an efficient recursive method published by Torben Andersen [1].

Andersen included psuedo-code in the original publication, this code presents a method providing improved speed and performance, as well as stability and some test cases.

Testing Method
The zernike values generated have been cross checked against hand calculated zernike value results for the first 4 radial order of Zernike polynomials, as well as against the results from the python package called Prysm.

The gradient values for the first 4 radial orders were similarly checked against hand calculated values and matched to 1E-8 accuracy. 

Languages Written
The algorithm has been written in python, and in rust, to provide a cross comparison of how to write the same algorithm in both languages.

References
[1] Torben B. Andersen, “Efficient and robust recurrence relations for the Zernike circle polynomials and their derivatives in Cartesian coordinates,” Opt. Express 26, 18878-18896 (2018)

Code
Note, for both the Rust and Python code, the ‘order’ or ‘index’ refers to the Ansi standard indexing, described here.

Both codes below are written under the MIT license, feel free to use however you would like, you do not need to reference me/this post for where you got the code (although I do hope you include the reference to the original paper!).

Usage
To use the Python code, simply navigate to the python directory, and run the 'recursive_zernike_generator.py'. The end of the python code has an example case built in that shows how to define the required inputs to the generator function, as well as how to call the function so it can be utilized in your own project.

To use the Rust code, navigate to the rust directory, and run 'cargo build' and then 'cargo run'. The code is provided in rust/src/main, where main() provides an example case for how to call and use the generator functions, and also a printout of the results. Calling cargo run will run the default test case and print out the results to your terminal. 
