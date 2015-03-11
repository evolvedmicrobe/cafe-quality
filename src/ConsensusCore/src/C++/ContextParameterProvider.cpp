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
    
    
    /* Autogenerated by OutputDiNucleotideContexts.py
     * Rows are Dark, Match, Stick (Branch is the reference)
     * Columns are Intercept, SNR, SNR^2, SNR^3 */
    Matrix<double> AA = {
        { 3.81720201515154, -0.400599603991922, 0.018614437316033, -0.000364751975465911 },
        { 3.1908276817261, 0.204011581063786, -0.0193345225294932, 0.000449003963971665 },
        { -1.25554885160522, 0.461682713184222, -0.040416513891308, 0.000914427334896904 }
    };
    Matrix<double> CC = {
        { -1.5430928904885, 1.61660832332086, -0.210997020100382, 0.00809223993524581 },
        { -3.74791833824717, 2.74983692727893, -0.313124643063006, 0.0113156595067565 },
        { -3.75220834760718, 1.75989172436986, -0.197306471137317, 0.00696163670699774 }
    };
    Matrix<double> GG = {
        { 8.81239057584291, -1.90170935230813, 0.182346682060559, -0.00539694591375577 },
        { 8.13157707315945, -1.15324215504543, 0.119778313278507, -0.00356129865508758 },
        { 7.25878582072985, -1.97560895808077, 0.185721972297463, -0.00527398059062249 }
    };
    Matrix<double> NA = {
        { 4.87271236251494, -1.0355068989864, 0.0610477738232239, -0.00130093910781944 },
        { 3.71545404320815, -0.111276259486349, 0.00298240855654832, -2.81445309541181e-05 },
        { 1.40041901429751, -0.435197508625283, 0.0255148449655911, -0.000581620461800338 }
    };
    Matrix<double> NC = {
        { 10.4302028387037, -3.27440991750123, 0.348001659986844, -0.0125520322736285 },
        { 7.92713496695142, -1.45779779769921, 0.161089563050961, -0.00592034192017757 },
        { 3.66876446637389, -1.11334654053442, 0.122005194567391, -0.00447662803608916 }
    };
    Matrix<double> NG = {
        { 2.02994606352135, -0.126177270318536, -0.0226746948052658, 0.00117260770884122 },
        { 3.12117495502696, 0.256982397810015, -0.0334213623483018, 0.00117900190125973 },
        { 0.443229415460136, 0.0313821995757453, -0.0222616337671001, 0.00102291356396239 }
    };
    Matrix<double> NT = {
        { 1.916265906613, -0.112957384644716, -0.00955897493038972, 0.000207132127635481 },
        { 2.44023498513114, 0.457060526378788, -0.0439016920522025, 0.00123404677896527 },
        { -2.08160499116939, 0.866536738681778, -0.106997300255227, 0.00380943870174409 }
    };
    Matrix<double> TT = {
        { 1.70503747459888, -0.0874041461144091, -0.012771794368405, 0.000528424406818803 },
        { 1.51229494096067, 0.513114778366272, -0.0569463299719687, 0.00190682385190301 },
        { -0.682649183406059, 0.153767976908084, -0.0421117527911171, 0.00193652817873698 }
    };
    
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
        sum = log(sum);
        
        double branch = -sum; // Branch probability is the reference, or 1 / sum
        
        // Now get the probabilities
        for(int i=0; i< 3; i++) {
            predicts[i] = log(predicts[i]) - sum;
        }
        TransitionParameters tp(predicts[1], predicts[2], branch, predicts[0]);
        return tp;
    }
    
}
