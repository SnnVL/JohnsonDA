# JohnsonDA

Python code accompanying the reference:

**Senne Van Loon and Steven J. Fletcher (2023). Foundations for Universal Non-Gaussian Data Assimilation, Geophysical Research Letters (accepted article)**

Please cite this work if you use this code.
@article{JohnsonVAR,
    title={{Foundations for Universal Non-Gaussian Data Assimilation}},
    author={Van Loon, Senne and Fletcher, Steven J.},
    journal={Geophysical Research Letters (accepted)},
    year={2023},
}


Files:

`run_L63.py` is the main run file, executing the 3DVAR algorithm with non-Gaussian errors

`plot_L63.py` plots the results from the run file

`mod_JohnsonDA.py` is the module file containing all necessary methods for the non-Gaussian 3DVAR

`models_rk.c` is the c-code integrating the Lorenz-63 model (should be compiled by running ```cc -fPIC -shared -o C_Lorenz.so models_rk4.c```)

Abstract:

*In many applications of data assimilation, especially when the size of the problem is large, a substantial assumption is made: all variables are well-described by Gaussian error statistics. This assumption has the advantage of making calculations considerably simpler, but it is often not valid, leading to biases in forecasts or, even worse, unphysical predictions. We propose a simple, but effective, way of replacing this assumption, by making use of transforming functions, while remaining consistent with Bayes' theorem. This method allows the errors to have any value of the skewness and kurtosis, and permits physical bounds for the variables. As such, the error distribution can conform better to the underlying statistics, reducing biases introduced by the Gaussian assumption. We apply this framework to a 3D variational data assimilation method, and find improved performance in a simple atmospheric toy model (Lorenz-63), compared to an all-Gaussian technique.*
