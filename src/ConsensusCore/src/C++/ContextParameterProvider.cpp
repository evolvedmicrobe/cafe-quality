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
        { -3.25348004644804, 0.109117295682203, -0.00478407086743064, 7.68229830247598e-05  },
        { -3.23879157079186, 0.121188914824026, -0.00715752186009457, 0.000160519218204073  },
        { -0.689416151193599, -0.466535535829279, 0.024542183986347, -0.000435458290662556  } };
    // Fit for context:  CC
    Matrix<double>  CC  = {
        { -5.20538947210907, 0.670790816932485, -0.0574130209298602, 0.00164629568878644  },
        { -1.58641491443679, -0.178410807934996, 0.0208028516863317, -0.000545250142713105  },
        { 0.940581728186671, -0.996712716219647, 0.0829954378705205, -0.00239928438398957  } };
    // Fit for context:  GG
    Matrix<double>  GG  = {
        { -3.1539834162978, 0.00227105478004632, 0.00375427292292107, -0.000129379476599567  },
        { -3.07063297042003, 0.162919555253349, -0.0168057857040007, 0.000507802181827188  },
        { -0.201010192793936, -0.577677765678952, 0.0340766005933718, -0.000656803578009733  } };
    // Fit for context:  TT
    Matrix<double>  TT  = {
        { -3.69995923007085, 0.168892331231276, -0.00230559944808306, -0.000217583682988837  },
        { -2.30708702574631, -0.00208655886285272, -0.00140804168812469, 6.2327169042168e-05  },
        { 0.676512170359799, -0.888008929374128, 0.0591267906332357, -0.00137707785866768  } };
    // Fit for context:  NA
    Matrix<double>  NA  = {
        { -2.7989451116997, 0.0800283597238306, -0.00379636553918088, 6.49417927985317e-05  },
        { -3.24147884207565, 0.0945318041196939, -0.00471742725753294, 9.38544188581566e-05  },
        { -0.146610214712297, -0.470967914226497, 0.0204419585309968, -0.000333016010149865  } };
    // Fit for context:  NC
    Matrix<double>  NC  = {
        { -4.91036493261879, 0.610796756954087, -0.0569008987583529, 0.00177463989325332  },
        { -2.48830865251312, 0.045683443720646, -0.00428870274457383, 0.000276929320038077  },
        { -2.58209716288575, -0.052148430028378, -0.0225639881328996, 0.00130123268129485  } };
    // Fit for context:  NG
    Matrix<double>  NG  = {
        { -3.45987744853226, 0.143922020016431, -0.00793991905740124, 0.000199286708772196  },
        { -2.68082863184399, 0.00461203484058528, 0.00140001163952265, -7.84441547496875e-05  },
        { -1.34079046453925, -0.291102281891434, 0.00639706321907262, 2.39270891924635e-05  } };
    // Fit for context:  NT
    Matrix<double>  NT  = {
        { -5.16586621525536, 0.645825449995914, -0.057169869144101, 0.00165876013435524  },
        { -2.62235602889202, -0.019022373597515, 0.00595499096486432, -0.000253124918958272  },
        { -1.31089221488273, -0.340095653242253, 0.0176648492740888, -0.00047397983472997  } };

    
    
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
