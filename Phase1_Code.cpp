#include <bits/stdc++.h>
using namespace std;


double objective_function1(double x);
double objective_function2(double x);
double objective_function3(double x);
double objective_function4(double x);
double objective_function5(double x);
double objective_function6(double x);

double numerical_Derivative(double f_xh0, double f_xh1, double h);
double numerical_SecondDerivative(double f_xh0, double f_x, double f_xh1, double h);

void Newton_Raphson(double a, double b, double (*objective_function)(double));
pair<double, double> Bounding_Phase(double a, double b, bool minimize, double (*objective_function)(double));

int c;
fstream out;

int main() {
    double a, b;       
    pair<double, double> result;
    int funcNum;

    cout << "Enter\n a = lower limit of x\n b = upper limit of x\n Enter 1 to minimize and 0 to maximize:\n";
    cin >> a >> b >> c;
    cout << "Enter Objective function number you want to use (1-6):\n";
    cin >> funcNum;

    /*Selecting Objective Function*/
    if (funcNum == 1) {
        /*Calling Bounding Phase Method*/
        result = (c == 1) ? Bounding_Phase(a, b, true, objective_function1) : Bounding_Phase(a, b, false, objective_function1);

        /*Calling Newton Raphson Method*/
        Newton_Raphson(result.first, result.second, objective_function1);
    }
    else if (funcNum == 2) {
        /*Calling Bounding Phase Method*/
        result = (c == 1) ? Bounding_Phase(a, b, true, objective_function2) : Bounding_Phase(a, b, false, objective_function2);

        /*Calling Newton Raphson Method*/
        Newton_Raphson(result.first, result.second, objective_function2);
    }
    else if (funcNum == 3) {
        /*Calling Bounding Phase Method*/
        result = (c == 1) ? Bounding_Phase(a, b, true, objective_function3) : Bounding_Phase(a, b, false, objective_function3);

        /*Calling Newton Raphson Method*/
        Newton_Raphson(result.first, result.second, objective_function3);
    }
    else if (funcNum == 4) {
        /*Calling Bounding Phase Method*/
        result = (c == 1) ? Bounding_Phase(a, b, true, objective_function4) : Bounding_Phase(a, b, false, objective_function4);

        /*Calling Newton Raphson Method*/
        Newton_Raphson(result.first, result.second, objective_function4);
    }
    else if (funcNum == 5) {
        /*Calling Bounding Phase Method*/
        result = (c == 1) ? Bounding_Phase(a, b, true, objective_function5) : Bounding_Phase(a, b, false, objective_function5);

        /*Calling Newton Raphson Method*/
        Newton_Raphson(result.first, result.second, objective_function5);
    }
    else if (funcNum == 6) {
        /*Calling Bounding Phase Method*/
        result = (c == 1) ? Bounding_Phase(a, b, true, objective_function6) : Bounding_Phase(a, b, false, objective_function6);

        /*Calling Newton Raphson Method*/
        Newton_Raphson(result.first, result.second, objective_function6);
    }
    else {
        cout << "Invalid Input\n";
    }

    return 0;
}


/*Objective Function*/
double objective_function1(double x) {
    return(pow(2*x-5, 4) - pow(x*x-1, 3));
}
double objective_function2(double x) {
    return(8 + pow(x, 3) - 2*x - 2*exp(x));
}
double objective_function3(double x) {
    return(4*x* sin(x));
}
double objective_function4(double x) {
    return(2*pow(x-3, 2) + exp(0.5*x*x));
}
double objective_function5(double x) {
    return(x*x - 10*exp(0.1*x));
}
double objective_function6(double x) {
    return(20*sin(x) - 15*x*x);
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

void Newton_Raphson(double a, double b, double (*objective_function)(double)) {
    out.open("Iterations.out", ios::app);
    out << ("\n**********************************\n");
	out << ("\nNewton Raphson Method\n");
    out << "#It\t\tx_(k)\t\tf(x_(k) - h)\t\tf(x_(k))\t\tf(x_(k) + h)\t\tf'(x_(k))\t\tf''(x_(k))\n";
    

    cout << ("\n**********************************\n");
	cout << ("Newton Raphson Method\n");
    

    /* Step 1: */
    uniform_real_distribution<double> unif(a, b);
    default_random_engine re;
    
    /* Getting a random double value between a and b*/
    double initial_guess = unif(re);   
    double epsilon = 0.0001;
    int k = 1, feval = 0;

    double x_k = initial_guess, x_k1, fk_Prime, fk_doublePrime;
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
    x_k1 = x_k - (fk_Prime / fk_doublePrime);

    f_xh0 = objective_function(x_k1 - h);
    f_xh1 = objective_function(x_k1 + h);
    feval = feval + 2;

    fk_Prime = numerical_Derivative(f_xh0, f_xh1, h);    


    /* Step 4: Termination condition*/
    if (fabs(fk_Prime) < epsilon) {
        cout << ((c == 1) ? "\nMinimum lies at " : "\nMaximum lies at ") << x_k1 << endl;
        cout << "Total number of function evaluations: " << feval << endl;

        /*Copying final results to file*/
        out << ("\n**********************************\n");
        out << ((c == 1) ? "\nMinimum lies at " : "\nMaximum lies at ") << x_k1 << endl;
        out << "Total number of function evaluations: " << feval << endl;
    }
    else {
        k = k + 1;
        x_k = x_k1;
        goto Step_2;    
    }
    out.close();
}

pair<double, double> Bounding_Phase (double a, double b, bool minimize, double (*objective_function)(double)) {
    out.open("Iterations.out", ios::out);
    out << ("\n**********************************\n");
	out << ("Bounding Phase Method\n");
    out << "#It\t\tx_(k-1)\t\tx_(k)\t\tx_(k+1)\t\tf(x_(k-1))\t\tf(x_(k))\t\tf(x_(k+1))\n";

    /* Step 1: */
    uniform_real_distribution<double> unif(a, b);
    default_random_engine re;
    
    double initial_guess = unif(re);  /* Getting a random double value between a and b*/
    double delta;
    int k = 0;
	
	cout << ("\n**********************************\n");
	cout << ("Bounding Phase Method\n");
	cout << ("Enter 'delta':\t");
    cin >> delta;


    /* Step 2: */
    double x1, x2, x3, f1, f2, f3;
    int feval;
    delta = fabs(delta);

    x2 = initial_guess; x1 = x2 - delta; x3 = x2 + delta; /*New points*/
	feval = 0; /*function evaluation*/
	f1 = objective_function(x1); /*Calculate objective_function*/
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
        if(f_k1 < f_k) {              /*If not terminated update the value of 'x' and f(x). Go to Step 3*/
            k = k + 1;
            x_k0 = x_k; x_k = x_k1;   
            f_k0 = f_k; f_k = f_k1;   
            goto Step_3;
        }
        else {
            cout << ("\n**********************************\n");
            cout << "Minimum lies between (" << x_k0 << ", " << x_k1 << ")" << endl;
            cout << "Total number of function evaluations: " << feval << endl;

            /*Copying final results to file*/
            out << ("\n**********************************\n");
            out << "Minimum lies between (" << x_k0 << ", " << x_k1 << ")" << endl;
            out << "Total number of function evaluations: " << feval << endl;
        }
    }

    /*Termination Condition for minimization problem*/
    else {
        if(f_k1 > f_k) {              /*If not terminated update the value of 'x' and f(x). Go to Step 3*/
            k = k + 1;
            x_k0 = x_k; x_k = x_k1;   
            f_k0 = f_k; f_k = f_k1;   
            goto Step_3;
        }
        else {
            cout << ("\n**********************************\n");
            cout << "Maximum lies between (" << x_k0 << ", " << x_k1 << ")" << endl;
            cout << "Total number of function evaluations: " << feval << endl;

            /*Copying final results to file*/
            out << ("\n**********************************\n");
            out << "Maximum lies between (" << x_k0 << ", " << x_k1 << ")" << endl;
            out << "Total number of function evaluations: " << feval << endl;
        }
    }  
    out.close();
    return make_pair(x_k0, x_k1);
}

