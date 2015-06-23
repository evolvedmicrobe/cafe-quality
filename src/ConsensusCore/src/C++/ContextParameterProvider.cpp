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
    
    
    Matrix<double>  AA  = {
        { -3.56366325147732, 0.213062470267352, -0.0147242069696605, 0.000351549429154497  },
        { -4.07477034464923, 0.259687926459921, -0.0130439403389001, 0.000197256427394713  },
        { 0.707509573270805, -0.828147896546462, 0.0524039113649921, -0.00105282871106833  } };
    // Fit for context:  CC
    Matrix<double>  CC  = {
        { -5.46541164550695, 0.731677480807461, -0.0665729662756421, 0.00199247921126783  },
        { -0.673194769919852, -0.55928583179955, 0.071612732415258, -0.00257195061120381  },
        { 1.69937115509728, -1.25910622753669, 0.109427539929014, -0.00307359016250297  } };
    // Fit for context:  GG
    Matrix<double>  GG  = {
        { -3.88637903232734, 0.231265622246946, -0.0186379412909804, 0.000417533559514484  },
        { -1.87221743170985, -0.227893294746469, 0.0195471474984407, -0.000487308703936223  },
        { 1.02192157753588, -1.02829324913349, 0.0847175848478058, -0.00230947739187885  } };
    // Fit for context:  TT
    Matrix<double>  TT  = {
        { -3.94050461752288, 0.245478039786449, -0.00986808751753838, 4.96609208843336e-05  },
        { -2.56085748845278, 0.0570397208303296, -0.00285173259458585, -2.28948684427145e-05  },
        { 2.31542340604299, -1.4922818622053, 0.126768962701795, -0.00364259350048365  } };
    // Fit for context:  NA
    Matrix<double>  NA  = {
        { -3.02575298026893, 0.135945701821264, -0.00871988369169161, 0.000188162974400923  },
        { -3.03478476634949, 0.0504400751870152, -0.00155501836394097, 2.91230346473706e-05  },
        { -0.285573322884753, -0.442159447043785, 0.0186741653746382, -0.000303086101398978  } };
    // Fit for context:  NC
    Matrix<double>  NC  = {
        { -6.18908075567119, 0.989594790506929, -0.0956503171462486, 0.00305959061984016  },
        { -1.88457739605028, -0.12275264979732, 0.0170359617053451, -0.000558883918986566  },
        { -2.10734932411511, -0.244286893265064, 0.00329029532422544, 0.000276740200101193  } };
    // Fit for context:  NG
    Matrix<double>  NG  = {
        { -3.44841280465889, 0.135543278883547, -0.0102370672949452, 0.000233236113268176  },
        { -2.33611759560326, -0.067325300456511, 0.00712595907402771, -0.000192028732220222  },
        { -1.08788438163626, -0.380095303545012, 0.0189786110583645, -0.000423479893703757  } };
    // Fit for context:  NT
    Matrix<double>  NT  = {
        { -3.96800324255753, 0.259069948477253, -0.0203032982450281, 0.000494731760604768  },
        { -3.02370710886002, 0.183089770632213, -0.0174583963432631, 0.000540431582190612  },
        { -1.51503259694717, -0.261982628866475, 0.010076855345182, -0.00022675901486311  } };
    
    
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
        
        double match = 1.0 / sum; // match probability is the reference, or 1 / sum
        //double branch = 1.0 / sum; // match probability is the reference, or 1 / sum
        
        // Now get the probabilities
        for(int i=0; i< 3; i++) {
            predicts[i] = predicts[i] / sum;
        }
        TransitionParameters tp(match, predicts[1], predicts[0], predicts[2]);
        //TransitionParameters tp(predicts[1], predicts[2], branch, predicts[0]);
        return tp;
    }
    
}
