//
//  ContextParameterProvider.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/27/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include "ContextParameterProvider.hpp"
using namespace std;


namespace ConsensusCore {
    
    SNR::SNR(double a, double c, double g, double t) : A(a), C(c), G(g), T(t) {}
    
    
    // SEE: /home/UNIXHOME/ndelaney/ccswork/Train/TrainLambdaModel.R
    // For how these parameters were generated.
    
    // Autogenerated by unitem package
    // Copy/Paste into ContextParameterProvider.cpp to test settings.
    // Fit for context:  AA
    Matrix<double>  AA  = {
        { 1.87273677034772, 1.87273677034772, -0.528748659234819, 0.0271376135739453, -0.000464386319247103  },
        { 1.786855033064, 1.78685503307358, -0.0231248695537108, -0.000287929406544749, 6.44630071493804e-06  },
        { 0.42707809960985, 0.427078099643392, -0.0243254424575994, -0.000233171208293363, -4.80615939419048e-05  } };
    // Fit for context:  CC
    Matrix<double>  CC  = {
        { 2.80434854719796, 2.80434854719796, -1.05280304678523, 0.0818471274882813, -0.00242565437497082  },
        { 2.02353123907654, 2.02353123907685, -0.0773581315055566, 0.0014766383410591, -0.000236642102138149  },
        { 1.58051867802742, 1.58051867802604, -0.697383792561278, 0.0712043801963417, -0.00253917030573885  } };
    // Fit for context:  GG
    Matrix<double>  GG  = {
        { 2.08218876147665, 2.08218876147665, -0.603047666731044, 0.0435301175644835, -0.00105732244064515  },
        { 1.79418977289188, 1.79418977289253, 0.0687586287533085, -0.00381939782675395, 8.36634321942822e-05  },
        { 1.15544512618762, 1.15544512617436, -0.501871658703042, 0.0412186766641648, -0.00105721602641759  } };
    // Fit for context:  TT
    Matrix<double>  TT  = {
        { 2.72511695915029, 2.72511695915029, -1.35060329442548, 0.11019312735876, -0.00324559724116098  },
        { 2.12960013125406, 2.12960013132381, -0.371080049227175, 0.0315172549358857, -0.000941263737682311  },
        { 1.19052060619789, 1.19052060619635, -0.608372283550254, 0.0420018278356971, -0.000998556224846127  } };
    // Fit for context:  NA
    Matrix<double>  NA  = {
        { 1.18378642815305, 1.18378642815305, -0.467223886779251, 0.0180038298001449, -0.000228366865213855  },
        { 1.60726372489734, 1.60726372489014, -0.0868354058082624, 0.00544065155236929, -0.000135271362827217  },
        { -0.0576847865984586, -0.0576847866068958, -0.0130873042912088, -0.00166393645054171, 1.35660733680042e-05  } };
    // Fit for context:  NC
    Matrix<double>  NC  = {
        { 2.98164569017579, 2.98164569017579, -1.73417821752682, 0.154708736003919, -0.00474415648726377  },
        { 1.94735356495702, 1.94735356499917, -0.183843549708076, 0.01818727256222, -0.000664710135761935  },
        { 1.21997341414346, 1.21997341414333, -0.671004413392223, 0.0705034497619124, -0.00249131126687544  } };
    // Fit for context:  NG
    Matrix<double>  NG  = {
        { 1.76454294928127, 1.76454294928127, -0.797044733183261, 0.047686662659488, -0.00106447289872615  },
        { 1.43054368767922, 1.43054368769032, 0.16199896475442, -0.0163886721381158, 0.000443325970021035  },
        { 0.119023809463059, 0.119023809268143, 0.0525469833016334, -0.0117698584896723, 0.000327755570412342  } };
    // Fit for context:  NT
    Matrix<double>  NT  = {
        { 2.6834950983038, 2.6834950983038, -1.45866147564357, 0.125966614629308, -0.00386791657567544  },
        { 1.70687055947819, 1.70687055945443, -0.0605565871910942, 0.0132012521351859, -0.000545418238333181  },
        { 0.664313535356815, 0.664313535248842, -0.233734495660121, 0.0199102448257229, -0.000664866857840956  } };

    
    static std::unordered_map<std::string, Matrix<double>* > parameter_store = { {"AA", &AA},{"CC", &CC},{"GG", &GG},{"NA", &NA},{"NC", &NC},{"NG", &NG},{"NT", &NT},{"TT", &TT}};
    
    TransitionParameters
    ContextParameterProvider::GetTransitionParameters(const string& context, const SNR& snrs)
    {
        auto params = *parameter_store[context];
        //Get the snr for the relevant channel
        auto channel  = context.at(1);
        double snr;
        switch(channel) {
                case 'A':
                    snr = snrs.A;
                    break;
                case 'C':
                    snr = snrs.C;
                    break;
                case 'G':
                    snr = snrs.G;
                    break;
                case 'T':
                    snr = snrs.T;
                    break;
            default:
                throw;
        }
        double snr2 = snr * snr;
        double snr3 = snr2 * snr;
        
        double predicts[3]; // Represents the XB portion
        double sum = 1.0;
        // Calculate each values contribution
        for(int i=0; i< 3; i ++) {
            auto xb = params[i][0] + snr * params[i][1] + snr2 * params[i][2] + snr3 * params[i][3];
            xb = exp(xb);
            predicts[i] = xb;
            sum += xb;
        }
        // Move to log space
        //sum = log(sum);
        
        double branch = 1.0 / sum; // Branch probability is the reference, or 1 / sum
        
        // Now get the probabilities
        for(int i=0; i< 3; i++) {
            predicts[i] = predicts[i] / sum;
        }
        TransitionParameters tp(predicts[1], predicts[2], branch, predicts[0]);
        return tp;
    }
    
}
