#include<bits/stdc++.h>

using namespace std;

double Actual_Function(const vector<double>& x);
double f(const vector<double>& x);
double g1(const vector<double>& x);
double g2(const vector<double>& x);
double g3(const vector<double>& x);
double g4(const vector<double>& x);
double g5(const vector<double>& x);
double g6(const vector<double>& x);
double bracket_Operator(const vector<double>& x, const function<double(const vector<double>&)>& func, double sigtau);
vector<double> Method_of_Multipliers (int n);

vector<double> gradient(const function<double(const vector<double>&)>& f, const vector<double>& x);
vector<vector<double>> inverseHessian(const function<double(const vector<double>&)>& f, const vector<double>& x);
vector<vector<double>> inverseMatrix(vector<vector<double>>& matrix, int n);

double norm(const vector<double>& x);
double dotProduct(const vector<double>& a, const vector<double>& b);
bool isLinearlyIndependent(const vector<double>& S, const vector<double>& grad);
vector<double> New_Direction(const vector<vector<double>>& H_inv, const vector<double>& grad);

double Objective_Function(double alpha);
double numerical_Derivative(double f_xh0, double f_xh1, double h);
double numerical_SecondDerivative(double f_xh0, double f_x, double f_xh1, double h);

pair<double, double> Bounding_Phase(double a, double b, bool minimize, double (*objective_function)(double));
double Newton_Raphson(double a, double b, double (*objective_function)(double));
vector<double> Newton_Multivariable_Gradient(double n);

vector<double> x, S;                    /* x is the current point, S is the search direction */
fstream out;                            /* Creating Output file */
int c = 1;                              /* For minimum c = 1, for maximum c = 0 */
int K, t;                                  /* Total number of iterations */
double final_feval = 0, R;              /* Total Function evalutions & penalty parameter */
vector<pair<double, double>> bounds;    /* Vector to store the bounds of each variable */
ifstream infile("input.txt");           /* Open Input file */

/* Vectors to store the greater than and equal to functions */
vector<function<double(const vector<double>&)>> greater_than_functions, equal_to_functions;  

/* Vectors to store the values of sigma and tau */
vector<double> sigma, tau;                                                                   


int main() {
    /* Read the number of variables from an input file */
    int n;
    infile >> n >> R;

    /* Read the bounds of each variable from an input file */
    for (int i = 0; i < n; i++) {
        double temp1, temp2;
        infile >> temp1 >> temp2;
        bounds.push_back(make_pair(temp1, temp2));
    }   

    /* Read the initial guess from an input file */
    for (int i = 0; i < n; ++i) {    
        double temp;
        infile >> temp;                 
        x.push_back(temp);
    }
    infile.close();

    /* Adding the greater than and equal to functions to the vector */
    /* Problem 1&2 */
    greater_than_functions.push_back(g1);
    greater_than_functions.push_back(g2);

    /* Problem 3 */
    // greater_than_functions.push_back(g3);   
    // greater_than_functions.push_back(g4);
    // greater_than_functions.push_back(g5);
    // greater_than_functions.push_back(g6);

    /* Calling Method of Multipliers */
    vector<double> solution = Method_of_Multipliers(n);


    /* Copying final results to file */
    out.open("Iterations.out", ios::app);
    out << ("\n********************************************\n");
    out << "Final Optimum Solution lies at:\n";
    for (auto i = 0; i < solution.size(); i++) {
        out << "x" << i+1 << " = " << setprecision(6) << solution[i] << "\n";
    }
    out << "\nValue at optimum point is: " << Actual_Function(solution);
    out << "\nTotal Function Evalution is: " << final_feval;
    out << "Total number of sequences: " << t << "\n";
    out << ("\n********************************************\n");
    out.close();

    /* Printing solution on the console */
    cout << ("\n********************************************\n");
	cout << ("Method of Multipliers (Constraint Optimisation Method)\n");
    cout << "Final Optimum Solution lies at:\n";
    for (auto i = 0; i < solution.size(); i++) {
        cout << "x" << i+1 << " = " << setprecision(6) << solution[i] << "\n";
    }
    cout << "\nValue at optimum point is: " << Actual_Function(solution) << endl;
    cout << "Total Function Evalution is: " << final_feval << "\n";
    cout << "Total number of sequences: " << t << "\n";
    cout << ("\n********************************************\n");
    return 0;
}


/* Actual objective function */
double Actual_Function(const vector<double>& x) {
    /* Problem 1 */
    return pow(x[0]-10, 3) + pow(x[1]-20, 3);

    /* Problem 2 */
    //return (-1 * pow(sin(2*3.14*x[0]), 3) * sin(2*3.14*x[1])) / (pow(x[0], 3)*(x[0]+x[1]));

    /* Problem 3 */
    // return x[0] + x[1] + x[2];
} 


/* Contraints */
double g1(const vector<double>& x) {
    /* Problem 1 */
    return pow(x[0]-5, 2) + pow(x[1]-5, 2) - 100;

    /* Problem 2 */
    //return x[1] - 1 - pow(x[0], 2);

    /* Problem 3 */
    // return (1 - 0.0025*(x[3]+x[5]));
}

double g2(const vector<double>& x) {
    /* Problem 1 */
    return 82.81 - pow(x[0]-6, 2) - pow(x[1]-5, 2);

    /* Problem 2 */
    //return x[0] - 1 - pow(x[1]-4, 2);

    /* Problem 3 */
    // return (1 - 0.0025*(x[4]+x[6]-x[3]));
}

double g3(const vector<double>& x) {
    /* Problem 3 */
    return (1 - 0.01*(x[7]-x[5]));
}

double g4(const vector<double>& x) {
    /* Problem 3 */
    return (x[0]*x[5]-833.33252*x[3]-100*x[0]+83333.333);
}

double g5(const vector<double>& x) {
    /* Problem 3 */
    return (x[1]*x[6]-1250*x[4]-x[1]*x[3]+1250*x[3]);
}

double g6(const vector<double>& x) {
    /* Problem 3 */
    return (x[2]*x[7]-1250000-x[2]*x[4]+2500*x[4]);
}


/* Bracket operator */
double bracket_Operator(const vector<double>& x, const function<double(const vector<double>&)>& func, double sigtau) {
    double sum = func(x) + sigtau;
    return (sum < 0)? sum : 0;
}


/* Penalty Function */
double f(const vector<double>& x) {
    double sum = 0;

    sum += Actual_Function(x);

    for (int i = 0; i < greater_than_functions.size(); ++i) {
        sum += R * (pow(bracket_Operator(x, greater_than_functions[i], sigma[i]), 2) - pow(sigma[i], 2));
    }
    for (int i = 0; i < equal_to_functions.size(); ++i) {
        sum += R * (pow(bracket_Operator(x, equal_to_functions[i], tau[i]), 2) - pow(tau[i], 2));
    }
    
    return sum;
}


/* Method of Multipliers Function */
vector<double> Method_of_Multipliers (int n) {
    
    /* Step 1: Initialisation */
    double epsilon1 = 1e-1, epsilon2 = 1e-1, prev_penalty = INT_MAX, curr_penalty = INT_MAX;
    vector<double> x_plus_1(n);
    t = 0;

    for (int i = 0; i < greater_than_functions.size(); ++i) {
        sigma.push_back(0);
    }
    for (int i = 0; i < equal_to_functions.size(); ++i) {
        tau.push_back(0);
    }

    prev_penalty = f(x);
    final_feval += 1;


    /* Step 2: Call Newton's Method Function*/
    Step_2:
    x_plus_1 = Newton_Multivariable_Gradient(n);
    

    /* Step 3: Termination Condition*/
    curr_penalty = f(x_plus_1);
    final_feval += 1;

    if (t != 0 && abs(prev_penalty - curr_penalty) <= epsilon2) {
        x = x_plus_1;
        return x;
    }

    /* Step 4: Update Sigma and Tau*/
    for (int i = 0; i < sigma.size(); ++i) {
        sigma[i] = bracket_Operator(x_plus_1, greater_than_functions[i], sigma[i]);
    }
    for (int i = 0; i < tau.size(); ++i) {
        tau[i] = bracket_Operator(x_plus_1, equal_to_functions[i], tau[i]);
    }

    prev_penalty = curr_penalty;
    x = x_plus_1;
    t++;
    goto Step_2;
}


/* Function to calculate the gradient of an objective function */
vector<double> gradient(const function<double(const vector<double>&)>& f, const vector<double>& x) {
    double h = 1e-4;                               // small constant for finite difference
    int n = x.size();
    vector<double> grad(n);

    for (int i = 0; i < n; ++i) {
        vector<double> x_plus_h = x;
        x_plus_h[i] += h;
        double f_plus_h = f(x_plus_h);

        vector<double> x_minus_h = x;
        x_minus_h[i] -= h;
        double f_minus_h = f(x_minus_h);

        grad[i] = (f_plus_h - f_minus_h) / (2 * h);
    }

    return grad;
}


/* Calculating the inverse Hessian matrix for a given objective function */
vector<vector<double>> inverseHessian(const function<double(const vector<double>&)>& f, const vector<double>& x) {
    double h = 1e-4;                            // small constant for finite difference
    int n = x.size();
    vector<vector<double>> H(n, vector<double>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            vector<double> x_plus_h = x;
            x_plus_h[i] += h;
            x_plus_h[j] += h;
            double f_plus_h = f(x_plus_h);

            vector<double> x_minus_h = x;
            x_minus_h[i] -= h;
            x_minus_h[j] -= h;
            double f_minus_h = f(x_minus_h);

            vector<double> x_i_plus_h_j_minus = x;
            x_i_plus_h_j_minus[i] += h;
            x_i_plus_h_j_minus[j] -= h;
            double f_i_plus_h_j_minus = f(x_i_plus_h_j_minus);

            vector<double> x_i_minus_h_j_plus = x;
            x_i_minus_h_j_plus[i] -= h;
            x_i_minus_h_j_plus[j] += h;
            double f_i_minus_h_j_plus = f(x_i_minus_h_j_plus);

            if(i != j) {
                H[i][j] = (f_plus_h - f_i_plus_h_j_minus - f_i_minus_h_j_plus + f_minus_h) / (4 * h * h);
            }
            else {
                H[i][j] = (f_plus_h - f_i_plus_h_j_minus - f_i_minus_h_j_plus + f_minus_h) / (h * h);
                H[j][i] = H[i][j];
            }
        }
    }

    /* Compute the inverse of the Hessian matrix */
    vector<vector<double>> H_inv = inverseMatrix(H, n);

    return H_inv;
}

/* Function to perform Gaussian elimination to find the inverse */
vector<vector<double>> inverseMatrix(vector<vector<double>>& matrix, int n) {
    vector<vector<double>> augmentedMatrix(n, vector<double>(2 * n, 0));

    /* Create an augmented matrix [A | I] where I is the identity matrix */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = matrix[i][j];
            if (i == j) augmentedMatrix[i][j + n] = 1;
        }
    }

    /* Perform Gaussian elimination */
    for (int i = 0; i < n; i++) {
        /* Partial pivoting to make the diagonal elements non-zero */
        int maxRow = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(augmentedMatrix[j][i]) > abs(augmentedMatrix[maxRow][i])) {
                maxRow = j;
            }
        }
        swap(augmentedMatrix[i], augmentedMatrix[maxRow]);

        double pivot = augmentedMatrix[i][i];
        if (pivot == 0.0) {
            cerr << "Matrix is singular. Inverse does not exist." << endl;
            exit(1);
        }

        for (int j = 0; j < 2 * n; j++) {
            augmentedMatrix[i][j] /= pivot;
        }

        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = augmentedMatrix[j][i];
                for (int k = 0; k < 2 * n; k++) {
                    augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                }
            }
        }
    }

    /* Extract the inverse matrix */
    vector<vector<double>> inverse(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmentedMatrix[i][j + n];
        }
    }

    return inverse;
}


/* Function to calculate the norm of a vector */
double norm(const vector<double>& x) {
    double sum = 0;
    for (double xi : x) {
        sum += xi * xi;
    }
    return sqrt(sum);
}


/* Function to calculate dot product of two vectors */
double dotProduct(const vector<double>& a, const vector<double>& b) {
    double sum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}


/* Function to check for linear independence between two directions */
bool isLinearlyIndependent(const vector<double>& S, const vector<double>& grad) {
    return (dotProduct(S, grad) <= 0)? true : false;
}


/* Function to calculate the new search direction */
vector<double> New_Direction(const vector<vector<double>>& H_inv, const vector<double>& grad) {
    vector<double> S;
    for (int i = 0; i < H_inv.size(); ++i) {
        double sum = 0;
        for (int j = 0; j < H_inv.size(); ++j) {
            sum += H_inv[i][j] * grad[j];
        }
        S.push_back(-sum);
    }
    return S;
}


/* Single variable function alpha */
double Objective_Function(double alpha) {
    vector<double> x_alpha(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        x_alpha[i] = x[i] - alpha * S[i];
    }
    return f(x_alpha);
}


/*First derivative of objective function*/
double numerical_Derivative(double f_xh0, double f_xh1, double h = 1e-6) {
    return (f_xh1 - f_xh0) / (2 * h);
}


/*Second derivative of objective function*/
double numerical_SecondDerivative(double f_xh0, double f_x, double f_xh1, double h = 1e-6) {
    return (f_xh1 - (2 * f_x) + f_xh0) / (h * h);
}


/*Newton Raphson Method*/
double Newton_Raphson(double a, double b, double (*objective_function)(double)) {
    out.open("Iterations.out", ios::app);
    out << ("\n**********************************\n");
	out << ("\nNewton Raphson Method\n");
    out << "#It\t\tx_(k)\t\tf(x_(k) - h)\t\tf(x_(k))\t\tf(x_(k) + h)\t\tf'(x_(k))\t\tf''(x_(k))\n";
    
    /* Step 1: */
    uniform_real_distribution<double> unif(a, b);
    default_random_engine re;
    
    /* Getting a random double value between a and b*/
    double initial_guess = unif(re);   
    double epsilon = 1e-1;
    int k = 1, feval = 0;

    double x_k = initial_guess, x_k1, fk_Prime, fk_doublePrime, prev_fk_Prime = 1e-6, diff = 1e-6;
    double f_x, f_xh0, f_xh1, h = 1e-6;

    f_xh0 = objective_function(x_k-h);
    f_xh1 = objective_function(x_k+h);
    f_x = objective_function(x_k);
    feval = feval + 3;

    fk_Prime = numerical_Derivative(f_xh0, f_xh1, h);


    /* Step 2: */
    Step_2:
    f_x = objective_function(x_k);
    feval++;
    fk_doublePrime = numerical_SecondDerivative(f_xh0, f_x, f_xh1, h);

    /*Copying values to file*/
    out << setprecision(5) << k << "\t\t" << setprecision(5) << x_k << "\t\t" << setprecision(5) << f_xh0 << "\t\t\t\t" << setprecision(5) << f_x << "\t\t\t" 
        << setprecision(5) << f_xh1 << "\t\t\t\t" << setprecision(5) << fk_Prime << "\t\t\t" << setprecision(5) << fk_doublePrime << "\n";


    /* Step 3: */
    if (fk_doublePrime == 0) {
        return x_k;
    }
    x_k1 = x_k - (fk_Prime / fk_doublePrime);

    f_xh0 = objective_function(x_k1 - h);
    f_xh1 = objective_function(x_k1 + h);
    feval = feval + 2;

    diff = fabs(fk_Prime - prev_fk_Prime);
    prev_fk_Prime = fk_Prime;
    fk_Prime = numerical_Derivative(f_xh0, f_xh1, h);    


    /* Step 4: Termination condition*/
if (fabs(fk_Prime) < epsilon || fabs(fk_Prime - prev_fk_Prime) == diff || fabs(fk_Prime - prev_fk_Prime) <= 1e-1 || fabs(x_k1 - x_k) <= 1e-1) {
        /*Copying final results to file*/
        out << ("\n**********************************\n");
        out << ((c == 1) ? "\nMinimum lies at " : "\nMaximum lies at ") << x_k1 << endl;
        out << "Total number of function evaluations: " << feval << endl;
        out.close();

        final_feval += feval;
        return x_k1;
    }
    else {
        k = k + 1;
        x_k = x_k1;
        goto Step_2;    
    }
    out.close();
}


/* Bounding Phase Method */
pair<double, double> Bounding_Phase (double a, double b, bool minimize, double (*objective_function)(double)) {
    out.open("Iterations.out", ios::app);
    out << ("\n**********************************\n");
	out << ("Bounding Phase Method\n");
    out << "#It\t\tx_(k-1)\t\tx_(k)\t\tx_(k+1)\t\tf(x_(k-1))\t\tf(x_(k))\t\tf(x_(k+1))\n";

    /* Step 1: */
    uniform_real_distribution<double> unif(a, b);
    default_random_engine re;
    
    double initial_guess = unif(re);             // Getting a random double value between a and b
    double delta = 0.1;
    int k = 0;
	

    /* Step 2: */
    double x1, x2, x3, f1, f2, f3;
    int feval;
    delta = fabs(delta);

    x2 = initial_guess; x1 = x2 - delta; x3 = x2 + delta;   // New points
	feval = 0;                                              // function evaluation
	f1 = objective_function(x1);                            // Calculate objective_function
	f2 = objective_function(x2);
	f3 = objective_function(x3);

	feval = feval + 3;

    if (minimize) {
        if (f1 >= f2 && f2 >= f3) {
            delta = fabs(delta);
        }
        else {
            delta = -fabs(delta);
        }
    }
    else {
        if (f1 >= f2 && f2 >= f3) {
            delta = -fabs(delta);
        }
        else {
            delta = fabs(delta);
        }
    }

    double x_k = x2, x_k1, x_k0 = x2; 
    double f_k = f2, f_k1, f_k0 = f2;


    /* Step 3: */
    Step_3:
	x_k1 = x_k + (pow(2, k) * delta);
    f_k1 = objective_function(x_k1);
    feval++;
    
    /*Copying values to file*/
    out << setprecision(5) << k << "\t\t" << setprecision(5) << x_k0 << "\t\t" << setprecision(5) << x_k << "\t\t" << setprecision(5) << x_k1 << "\t\t" 
        << setprecision(5) << f_k0 << "\t\t\t" << setprecision(5) << f_k << "\t\t\t" << setprecision(5) << f_k1 << "\n";


    /* Step 4: */

    /*Termination Condition for minimization problem*/
    if(minimize) {
        if(f_k1 < f_k) {                // If not terminated update the value of 'x' and f(x). Go to Step 3
            k = k + 1;
            x_k0 = x_k; x_k = x_k1;   
            f_k0 = f_k; f_k = f_k1;   
            goto Step_3;
        }
        else {                          // If terminated, print the final results to file
            out << ("\n**********************************\n");
            out << "Minimum lies between (" << x_k0 << ", " << x_k1 << ")" << endl;
            out << "Total number of function evaluations: " << feval << endl;

            final_feval += feval; 
        }
    }

    /*Termination Condition for minimization problem*/
    else {
        if(f_k1 > f_k) {                // If not terminated update the value of 'x' and f(x). Go to Step 3
            k = k + 1;
            x_k0 = x_k; x_k = x_k1;   
            f_k0 = f_k; f_k = f_k1;   
            goto Step_3;
        }
        else {                          // If terminated, print the final results to file
            out << ("\n**********************************\n");
            out << "Maximum lies between (" << x_k0 << ", " << x_k1 << ")" << endl;
            out << "Total number of function evaluations: " << feval << endl;

            final_feval += feval;
        }
    }  
    out.close();
    return make_pair(x_k0, x_k1);
}


/* Newton's method - Gradient based Multivariable optimization algorithm */
vector<double> Newton_Multivariable_Gradient(double n) {

    /* Step 1: Initialisation */
    double M = 10, epsilon1 = 1e-1, epsilon2 = 1e-1;
    K = 0;

    vector<double> x_plus_1(n);

    /* Copying initial guess to file */
    out.open("Iterations.out", ios::app);
    out << ("\n**********************************\n");
	out << ("\nNewton Method (Gradient based Multivariable optimization Algorithm)\n");
    out << "#Iterations\tCurrent Iteration points\n";
    for (double xi : x) {
        out << setprecision(6) << xi << "\t";
    }
    out << "\n";
    out.close();


    /* Step 2: */
    Step_2:
    vector<double> grad = gradient(f, x);
    final_feval += 2 * n;

    /* Step 3: */
    Step_3:
    double norm_grad = norm(grad);

    /* Termination condition 1: Based on gradient and iterations */
    if (norm_grad < epsilon1 || K > M) {
        return x;
    }


    /* Step 4: */
    vector<vector<double>> H_inv = inverseHessian(f, x);
    final_feval += (2 * n * n) + 1; 
    S = New_Direction(H_inv, grad);
    
    if (isLinearlyIndependent(S, grad)) {            // checking for linear independence
  
        for (auto i = 0; i < n; ++i) {
            pair<double, double> temp = Bounding_Phase(bounds[i].first, bounds[i].second, true, Objective_Function);
            double alpha = Newton_Raphson(temp.first, temp.second, Objective_Function);
            x_plus_1[i] = (x[i] - alpha * S[i]);
        }

        /* Copying values to file */
        out.open("Iterations.out", ios::app);
        out << "\n";   
        for (double xi : x_plus_1) {
            out << setprecision(6) << xi << "\t";
        }   
        out << "\n";
        out.close();
    }
    else {
        /* Restarting the algorithm */
        for (auto i = 0; i < n; ++i) {
            x[i] = 16;
        }
        goto Step_2;
    }

    /* Termination condition 2 : Based on unidirectional search */
    vector<double> grad_plus_1 = gradient(f, x_plus_1);
    final_feval += 2 * n;

    if (abs(dotProduct(grad_plus_1, grad)) < epsilon2) {
        return x_plus_1;
    }

    /* Step 5: */
    vector<double> x_New_minus_Old (x.size());
    for (size_t i = 0; i < n; ++i) {
        x_New_minus_Old[i] = x_plus_1[i] - x[i];
    }

    /* Termination condition 3 : Based on norm of current and previous iteration points */
    if ((norm(x_New_minus_Old)/norm(x)) < epsilon1) {
        return x_plus_1;
    }
    else {
        x.clear();
        x = x_plus_1;
        grad = grad_plus_1;
        K++;
        goto Step_3;
    }
}