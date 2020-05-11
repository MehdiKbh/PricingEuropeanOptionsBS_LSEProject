#include "MonteCarloEngine.h"
#include "GBM.h"
#include "Payoffs.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

//Constructor
MC_engine_GBM::MC_engine_GBM(GBM S, double Maturity, int nSimuls){
    model = S;
    T = Maturity;
    nbSimuls = nSimuls;

}

//Classical MC estimator
double MC_engine_GBM::computeEuropeanPrice( double (*Payoff)(double spot, double* params), double* params_payoff, double* weights )
{

    double price = 0;
    int n = this->nbSimuls;
    double* simuls = new double [n]; //Keep the simulations in memory  in order to compute the sample variance
    Matrix w(weights, model.get_N()); //Weights in a column vector

    for(int i=0; i<n; i++){
        Matrix simul_i = model.samplepath(T, 1).get_exctract(0, this->model.get_N()-1, 1, 1); //Simulation of samples of S(T) (column vector)
        double X = (w.T()*simul_i).get_element(0,0); //Weighted average of w*S (basket of stocks)
        simuls[i] = exp(-model.get_r()*T)*(*Payoff)(X, params_payoff); //Compute the payoff for this simulation and store it in memory
        price += simuls[i]; //Update the mean
    }
    price /= (float)n; //Average

    //Compute sample variance using the simulations stored in memory
    double sample_var=0;
    for(int i=0; i<n; i++){
        sample_var += (simuls[i] - price)*(simuls[i] - price);
    }
    sample_var /= (float)(n-1);

    //Print the confidence interval
    cout << "95% confidence interval of the option's price:\n";
    cout << "(" << price - 1.96*(sqrt(sample_var/(float)n)) << "," << price + 1.96*(sqrt(sample_var/(float)n)) << ")\n";

    delete simuls;
    return price;
}



//Control variate estimator
double MC_engine_GBM::computeEuropeanPrice_controlVariate( double (*Payoff)(double spot, double* params), double* params_payoff, double* weights )
{

    double mean_Y = 0;
    double mean_X = 0;
    int n = this->nbSimuls;
    double* simuls_Y = new double [n]; //Keep simulations in memory for computing sample variance
    double* simuls_X = new double [n];
    Matrix w(weights, model.get_N());

    for(int i=0; i<n; i++){
        Matrix simul_i = model.samplepath(T, 1).get_exctract(0, this->model.get_N()-1, 1, 1); //Simulation of samples of S(T) (column vector)
        double X = (w.T()*simul_i).get_element(0,0);//Weighted average of w*S (Basket of stocks)


        simuls_X[i] = exp(-model.get_r()*T)*X; //Compute the discounted value of the basket and store it in memory
        simuls_Y[i] = exp(-model.get_r()*T)*(*Payoff)(X, params_payoff); //Compute the payoff for this simulation and store it in memory

        //Update the means
        mean_X += simuls_X[i];
        mean_Y += simuls_Y[i];
    }
    //Average
    mean_X /= (float)n;
    mean_Y /= (float)n;

    //Compute b
    double var_X = 0;
    double cov_XY = 0;
    for(int j=0; j<n; j++){
        var_X += (simuls_X[j] - mean_X)*(simuls_X[j] - mean_X);
        cov_XY += (simuls_X[j] - mean_X)*(simuls_Y[j] - mean_Y);
    }
    double b = cov_XY/var_X;

    //Compute the value of the estimate
    double X0 = (w.T()*this->model.get_S0()).get_element(0,0);
    double price = mean_Y - b*(mean_X - X0);

    //Compute sample variance using the simulations stored in memory
    double sample_var=0;
    for(int i=0; i<n; i++){
        sample_var += ( (simuls_Y[i] - b*(simuls_X[i] - X0)) - price)*((simuls_Y[i] - b*(simuls_X[i] - X0)) - price);
    }
    sample_var /= (float)(n-1);

    //Print the confidence interval
    cout << "95% confidence interval of the option's price:\n";
    cout << "(" << price - 1.96*(sqrt(sample_var/(float)n)) << "," << price + 1.96*(sqrt(sample_var/(float)n)) << ")\n";

    delete simuls_X;
    delete simuls_Y;
    return price;
    }


//Test prices for call and put options for the classical MC estimator
void MC_engine_GBM::testCallPut(double* K, double* weights, double* res ){

    //Call option
    double price_call = 0;
    int n = this->nbSimuls;
    cout << n << ",";
    double* simuls_call = new double [n];
    Matrix w(weights, model.get_N());

    for(int i=0; i<n; i++){
        //Simulation of samples of S(T) (vector)
        Matrix simul_i_call = model.samplepath(T, 1).get_exctract(0, this->model.get_N()-1, 1, 1);
        //Weighted average of w*S
        double X = (w.T()*simul_i_call).get_element(0,0);
        //Compute the payoff for this simulation and store it in memory
        simuls_call[i] = exp(-model.get_r()*T)*callPayoff(X, K);
        //Update the mean
        price_call += simuls_call[i];
    }
    price_call /= (float)n;
    cout << price_call << ",";

    //Compute sample variance for the call
    double sample_var_call=0;
    for(int i=0; i<n; i++){
        sample_var_call += (simuls_call[i] - price_call)*(simuls_call[i] - price_call);
    }
    sample_var_call /= (float)(n-1);

    double call_low = price_call - 1.96*(sqrt(sample_var_call/(float)n));
    double call_up = price_call + 1.96*(sqrt(sample_var_call/(float)n));
    cout << call_low << "," << call_up << ",";


    //Put option
    double price_put = 0;
    double* simuls_put = new double [n];

    for(int i=0; i<n; i++){
        //Simulation of samples of S(T) (vector)
        Matrix simul_i_put = model.samplepath(T, 1).get_exctract(0, this->model.get_N()-1, 1, 1);
        //Weighted average of w*S
        double X = (w.T()*simul_i_put).get_element(0,0);
        //Compute the payoff for this simulation and store it in memory
        simuls_put[i] = exp(-model.get_r()*T)*putPayoff(X, K);
        //Update the mean
        price_put += simuls_put[i];
    }
    price_put /= (float)n;

     //Compute sample variance for the put
    double sample_var_put=0;
    for(int i=0; i<n; i++){
        sample_var_put += (simuls_put[i] - price_put)*(simuls_put[i] - price_put);
    }
    sample_var_put /= (float)(n-1);

    double put_low = price_put - 1.96*(sqrt(sample_var_put/(float)n));
    double put_up = price_put + 1.96*(sqrt(sample_var_put/(float)n));
    cout << put_low << "," << put_up << ",";


    //C-P
    cout << price_call - price_put << ",";
    double sample_var_cp=0;
    for(int i=0; i<n; i++){
        sample_var_cp += ((simuls_call[i] - simuls_put[i]) - (price_call - price_put))*((simuls_call[i] - simuls_put[i]) - (price_call - price_put));
    }
    sample_var_cp /= (float)(n-1);


    double cp_low = (price_call-price_put) - 1.96*(sqrt(sample_var_cp/(float)this->nbSimuls));
    double cp_up = (price_call-price_put)+ 1.96*(sqrt(sample_var_cp/(float)this->nbSimuls));
    cout << cp_low << "," << cp_up << endl;

    delete simuls_call;
    delete simuls_put;

}



//Test prices for call and put options for the control variate estimator
void MC_engine_GBM::testCallPut_control(double *K, double* weights){


    //Call price
    double mean_Y_call = 0;
    double mean_X_call = 0;
    int n = this->nbSimuls;
    cout << n << ",";
    double* simuls_Y_call = new double [n];
    double* simuls_X_call = new double [n];
    Matrix w(weights, model.get_N());

    for(int i=0; i<n; i++){
        //Simulation of samples of S(T) (column vector)
        Matrix simul_i = model.samplepath(T, 1).get_exctract(0, this->model.get_N()-1, 1, 1);
        //Weighted average of w*S
        double X = (w.T()*simul_i).get_element(0,0);
        //Compute discounted basket price
        simuls_X_call[i] = exp(-model.get_r()*T)*X;
        //Compute the payoff for this simulation and store it in memory
        simuls_Y_call[i] = exp(-model.get_r()*T)*callPayoff(X, K);
        //Update the mean
        mean_X_call += simuls_X_call[i];
        mean_Y_call += simuls_Y_call[i];
    }
    mean_X_call /= (float)n;
    mean_Y_call /= (float)n;

    //Compute b for the call option
    double var_X_call = 0;
    double cov_XY_call = 0;
    double var_Y_call = 0;
    for(int j=0; j<n; j++){
        var_X_call += (simuls_X_call[j] - mean_X_call)*(simuls_X_call[j] - mean_X_call);
        cov_XY_call += (simuls_X_call[j] - mean_X_call)*(simuls_Y_call[j] - mean_Y_call);
        var_Y_call += (simuls_Y_call[j] - mean_Y_call)*(simuls_Y_call[j] - mean_Y_call);
    }
    double b_call = cov_XY_call/var_X_call;

    double X0 = (w.T()*this->model.get_S0()).get_element(0,0);
    double price_call = mean_Y_call - b_call*(mean_X_call - X0);

    //Compute sample variance using the simulations stored in memory
    double sample_var_call=0;
    for(int i=0; i<n; i++){
        sample_var_call += ( (simuls_Y_call[i] - b_call*(simuls_X_call[i] - X0)) - price_call)*((simuls_Y_call[i] - b_call*(simuls_X_call[i] - X0)) - price_call);
    }
    sample_var_call /= (float)(n-1);
    double down_call = price_call - 1.96*(sqrt(sample_var_call/(float)n));
    double up_call =  price_call + 1.96*(sqrt(sample_var_call/(float)n));
    cout << price_call << "," << down_call << "," << up_call << ",";


    //Put price
    double mean_Y_put = 0;
    double mean_X_put = 0;
    double* simuls_Y_put = new double [n];
    double* simuls_X_put = new double [n];

    for(int i=0; i<n; i++){
        //Simulation of samples of S(T) (vector)
        Matrix simul_i = model.samplepath(T, 1).get_exctract(0, this->model.get_N()-1, 1, 1);
        //Weighted average of w*S
        double X = (w.T()*simul_i).get_element(0,0);
        //Discounted basket price
        simuls_X_put[i] = exp(-model.get_r()*T)*X;
        //Compute the payoff for this simulation and store it in memory
        simuls_Y_put[i] = exp(-model.get_r()*T)*putPayoff(X, K);
        //Update the mean
        mean_X_put += simuls_X_put[i];
        mean_Y_put += simuls_Y_put[i];
    }
    mean_X_put /= (float)n;
    mean_Y_put /= (float)n;

    //Compute b for the put price
    double var_X_put = 0;
    double cov_XY_put = 0;
    double var_Y_put = 0;
    for(int j=0; j<n; j++){
        var_X_put += (simuls_X_put[j] - mean_X_put)*(simuls_X_put[j] - mean_X_put);
        cov_XY_put += (simuls_X_put[j] - mean_X_put)*(simuls_Y_put[j] - mean_Y_put);
        var_Y_put += (simuls_Y_put[j] - mean_Y_put)*(simuls_Y_put[j] - mean_Y_put);
    }
    double b_put = cov_XY_put/var_X_put;

    double price_put = mean_Y_put - b_put*(mean_X_put - X0);

    //Compute sample variance using the simulations stored in memory
    double sample_var_put=0;
    for(int i=0; i<n; i++){
        sample_var_put += ( (simuls_Y_put[i] - b_put*(simuls_X_put[i] - X0)) - price_put)*((simuls_Y_put[i] - b_put*(simuls_X_put[i] - X0)) - price_put);
    }
    sample_var_put /= (float)(n-1);
    double down_put = price_put - 1.96*(sqrt(sample_var_put/(float)n));
    double up_put =  price_put + 1.96*(sqrt(sample_var_put/(float)n));
    cout << price_put << "," << down_put << "," << up_put << ",";


    //C-P
    double c_p = price_call - price_put;
    cout << c_p << ",";
    double var_cp = 0;
    for(int i=0; i<n; i++){
        double call_i = simuls_Y_call[i] - b_call*(simuls_X_call[i] - X0);
        double put_i = simuls_Y_put[i] - b_put*(simuls_X_put[i] - X0);
        double cp_i = call_i - put_i;

        var_cp += ( cp_i - c_p )*( cp_i - c_p );
    }
    var_cp /= (float)(n-1);

    double down_cp = (price_call - price_put) - 1.96*(sqrt(var_cp/(float)n));
    double up_cp = (price_call - price_put) + 1.96*(sqrt(var_cp/(float)n));
    cout << down_cp << "," << up_cp << "," << cov_XY_call/sqrt(var_X_call*var_Y_call) << "," << cov_XY_put/sqrt(var_X_put*var_Y_put) << endl;

}


