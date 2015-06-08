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
        { -3.57741855801257, 0.00284057592827521, 0.00487752334278775, -0.000164210962674474  },
        { -3.09244317485133, 0.0368499858770809, -0.00283838380254529, 5.9242137281507e-05  },
        { 0.149420291744396, -0.565502169337908, 0.0325218874456029, -0.000593081138032452  } };
    // Fit for context:  CC
    Matrix<double>  CC  = {
        { -4.94272678453286, 0.549132261343589, -0.0547112795170815, 0.00170518697330163  },
        { -1.2707282376388, -0.621298951028855, 0.0801305916566298, -0.00301479171083435  },
        { 1.09504682733958, -0.975541103137689, 0.0795201493284305, -0.00203407738696825  } };
    // Fit for context:  GG
    Matrix<double>  GG  = {
        { -3.87377311027209, 0.0902772025622097, -0.00396275248815696, -2.77908476778987e-05  },
        { -2.6780977927035, -0.094707198355993, 0.00524981352158013, -0.000186986375017649  },
        { 0.520618032890895, -0.803436236559507, 0.0610740391697061, -0.00155013480663096  } };
    // Fit for context:  TT
    Matrix<double>  TT  = {
        { -3.76626093932139, 0.15652826314597, -0.00651912926466681, 6.49242036189152e-05  },
        { -3.41057727247402, 0.250202519308169, -0.0355597646872119, 0.00122887770498529  },
        { 0.832717865382825, -0.886109040849563, 0.0632838037475672, -0.00148984965053754  } };
    // Fit for context:  NA
    Matrix<double>  NA  = {
        { -3.25696710131804, 0.0456286158584122, -0.00163110998806611, 3.94576992886561e-05  },
        { -3.4331078149566, 0.0263634723284365, 0.000712256880125147, -5.55624760288029e-05  },
        { -0.455329686650815, -0.373456892695384, 0.0133560078384758, -0.000130728038623462  } };
    // Fit for context:  NC
    Matrix<double>  NC  = {
        { -4.47019697584345, 0.275689971004247, -0.017390576657798, 0.000291812514785925  },
        { -2.95462818681068, 0.0366753950312157, -0.000138550762629993, -0.000120084150910541  },
        { -2.14552806543074, -0.215508722091391, 0.0103713863138445, -0.000269155986232906  } };
    // Fit for context:  NG
    Matrix<double>  NG  = {
        { -3.33047198391906, 0.00982242878996712, 8.25782971326983e-05, -3.02502254908397e-07  },
        { -2.69368130908551, -0.0922582515445852, 0.00439969713856122, -4.46086544327924e-05  },
        { -0.924986427325525, -0.396196190930418, 0.0212632882434381, -0.000423676445459592  } };
    // Fit for context:  NT
    Matrix<double>  NT  = {
        { -3.17731160857182, -0.0426161699551311, 0.00607843567208108, -0.000208835221989087  },
        { -2.69891453107569, -0.0906019426740408, 0.00570884560011702, -8.70025818462857e-05  },
        { -1.5622974475084, -0.206339192033753, 0.00564379522803271, -5.08674424310348e-05  } };
    
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
