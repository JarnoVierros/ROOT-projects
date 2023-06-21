
string create_interval(float start, float stop) {
    return "((("+to_string(start)+"<x) ? 1 : 0) - (("+to_string(stop)+"<x) ? 1 : 0))";
}

string get_par(TF1 fit, int par) {
    return to_string(fit.GetParameter(par));
}

class Particle {
    public:
        //pion=1, rho=2
        int type;
        float p;
        float p_t;
        int charge;
        float eta;
        float phi;
        float p_x;
        float p_y;
        float p_z;
        float m;
        float E;

        Particle(int type_in) {
            type = type_in;   
        }

        void calculate_3d_momentum() {
            if (p_t==0||eta==0||phi==0) {
                throw invalid_argument("required variables not set");
            }
            p_x = p_t*cos(phi);
            p_y = p_t*sin(phi);
            p_z = p_t*sinh(eta);
        }

        void calculate_energy() {
            E = sqrt(pow(p, 2) + pow(m, 2));
        }

        void calculate_total_momentum() {
            p = sqrt(pow(p_x, 2)+pow(p_y, 2)+pow(p_z, 2));
        }

        void calculate_mass() {
            m = sqrt(pow(E, 2) - pow(p, 2));
        }
};

class Event {
    public:
        Particle** particles;
        int particle_count;
        float ThxR;
        float ThyR;
        float ThxL;
        float ThyL;

        Event(Particle* particles_in[], int particle_count_in) {
            particles = particles_in;
            particle_count = particle_count_in;
        }

        int calculate_total_charge() {
            float total_charge = 0.;
            for (int i=0; i<particle_count; ++i) {
                Particle* particle = *(particles+i);
                total_charge += particle->charge;
            }
            return total_charge;
        }

        void reconstruct_2_rhos(Particle* rhos[2][2]) {
            if (particle_count != 4) {
                throw invalid_argument("invalid number of particles");
            }

            Particle* positives[2];
            int positive_index = 0;
            Particle* negatives[2];
            int negative_index = 0;

            for (int i=0; i<4; ++i) {
                Particle* current_particle = *(particles+i);
                if (current_particle->charge == 1) {
                    positives[positive_index] = current_particle;
                    positive_index++;
                } else if (current_particle->charge == -1) {
                    negatives[negative_index] = current_particle;
                    negative_index++;
                } else {
                    throw invalid_argument("particle has invalid charge");
                }
            }

            Particle* rho1_c1 = reconstruct_rho(positives[0], negatives[0]);
            rhos[0][0] = rho1_c1;
            Particle* rho2_c1 = reconstruct_rho(positives[1], negatives[1]);
            rhos[0][1] = rho2_c1;
            
            Particle* rho1_c2 = reconstruct_rho(positives[0], negatives[1]);
            rhos[1][0] = rho1_c2;
            Particle* rho2_c2 = reconstruct_rho(positives[1], negatives[0]);
            rhos[1][1] = rho2_c2;
        }
    
        Particle* reconstruct_rho(Particle* positive, Particle* negative) {
            Particle* rho = new Particle(2);
            rho->E = positive->E + negative->E;
            rho->p_x = positive->p_x + negative->p_x;
            rho->p_y = positive->p_y + negative->p_y;
            rho->p_z = positive->p_z + negative->p_z;
            rho->calculate_total_momentum();
            rho->calculate_mass();
            return rho;
        }
};


void ntuple_analysis_background_fit() {

    const string filename = "TOTEM43.root"; //"TOTEM43.root", kpkm.roo, 110000.root
    const float muon_mass = 139.57039;
    const float rho_mass = 730; //770

    //gROOT->SetStyle("Plain");   // set plain TStyle
    gStyle->SetOptStat(0); // draw statistics on plots,
                                // (0) for no output
    //gStyle->SetOptFit(1111);    // draw fit results on plot,
                                // (0) for no ouput
    //gStyle->SetPalette(57);     // set color map
    //gStyle->SetOptTitle(0);       // suppress title box

    TFile *file = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)file->Get("tree");
    TTreeReader Reader(tree);
    //TTreeReader Reader("tree", file);

    auto rho_masses = new TH2F("rho_masses", "Masses of potential rho particles;m1/MeV;m2/MeV",200,200,1400,200,200,1400);

    TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
    TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");
    TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
    TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
    TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");

    TTreeReaderValue<Float_t> ThxR(Reader, "ThxR");
    TTreeReaderValue<Float_t> ThyR(Reader, "ThyR");
    TTreeReaderValue<Float_t> ThxL(Reader, "ThxL");
    TTreeReaderValue<Float_t> ThyL(Reader, "ThyL");

    while (Reader.Next()) {

        if (trk_p.GetSize() != 4) {
            continue;
        }

        Particle* particles[4];

        for (int i=0;i<4;++i) {
            Particle* particle = new Particle(1);
            particle->p = 1000*trk_p[i];
            particle->p_t = 1000*trk_pt[i];
            particle->charge = trk_q[i];
            particle->eta = trk_eta[i];
            particle->phi = trk_phi[i];
            particle->m = muon_mass;
            particle->calculate_3d_momentum();
            particle->calculate_energy();
            particles[i] = particle;
        }
        
        /*
        for (Particle* particle : particles) {
            //cout << "momentum: " << particle.p_t << ", charge: " << particle.charge << endl;
            cout << "p_x: " << particle->p_x << ", p_y: " << particle->p_y << ", p_z: " << particle->p_z << ", Dp: " 
            << particle->p - sqrt(pow(particle->p_x, 2)+pow(particle->p_y, 2)+pow(particle->p_z, 2)) 
            << ", Dp_t: " << particle->p_t - sqrt(pow(particle->p_x, 2)+pow(particle->p_y, 2)) << ", charge: " << particle->charge << endl;
        }
        */

        Event current_event(particles, 4);

        current_event.ThxR = *ThxR;
        current_event.ThyR = *ThyR;
        current_event.ThxL = *ThxL;
        current_event.ThyL = *ThyL;

        float total_charge = current_event.calculate_total_charge();
        if (total_charge != 0) {
            //cout << "INVALID" << endl << endl;
            continue;
        }

        Particle* rhos[2][2];
        current_event.reconstruct_2_rhos(rhos);

        float masses[2];

        for (int i=0; i<2; ++i) {
            for (int j=0; j<2; ++j) {
                Particle* rho = rhos[i][j];
                masses[j] = rho->m;
                //cout << "m: " << rho->m << endl;
            }
            rho_masses->Fill(masses[0], masses[1]);
        }

        //cout << endl;
    }


    auto main = new TCanvas("Canvas1","Canvas1");
    rho_masses->Draw("Colz");

    TLine line1 = TLine(rho_mass, 200, rho_mass, 1400);
    line1.DrawClone();

    TLine line2 = TLine(200, rho_mass, 1400, rho_mass);
    line2.DrawClone();


    TArrow arrow1(850,550,750,710,0.02,"|>");
    arrow1.SetLineWidth(3);
    arrow1.DrawClone();

    TLatex text1(860,500,"\\rho");
    text1.DrawClone();


    TArrow arrow2(350,350,450,450,0.02,"|>");
    arrow2.SetLineWidth(3);
    arrow2.DrawClone();

    TLatex text2(290,300,"K^{0}_{s}");
    text2.DrawClone();


    main->Print("rho_masses.pdf","RECREATE");

    auto projections = new TCanvas("Canvas2","Canvas2");
    projections->Divide(1,2);

    projections->cd(1);
    rho_masses->ProjectionX()->SetTitle("Mass of first potential rho particle");
    rho_masses->ProjectionX()->GetYaxis()->SetTitle("events");
    rho_masses->ProjectionX()->Draw();
    //TLine line3 = TLine(rho_mass, 0, rho_mass, 1.05*rho_masses->ProjectionX()->GetMaximum());
    //line3.DrawClone();

    const float rho_fit_min = 629.; //629.535 //650
    const float rho_fit_max = 850.;
    const float k_fit_min = 470.;
    const float k_fit_max = 530.;

    TString linear_bckg_string = "("+create_interval(k_fit_max, rho_fit_min)+"+"+create_interval(0., k_fit_min)+")*([0]+[1]*x)";
    cout << linear_bckg_string << endl;
    TF1 linear_bckg("linear_bckg", linear_bckg_string, 200, 1400);
    linear_bckg.SetParameters(1, 1);

    rho_masses->ProjectionX()->Fit(&linear_bckg, "","",400,629.535);


    TString k_fit_string = get_par(linear_bckg, 0)+"+"+get_par(linear_bckg, 1)+"*x+[0]*TMath::Gaus(x,[1],[2])";
    cout << k_fit_string << endl;
    TF1 k_fit("k_fit", k_fit_string, 200, 1400);
    k_fit.SetParameters(8.19813e+02, 4.97588e+02, 1.19288e+01);

    rho_masses->ProjectionX()->Fit(&k_fit, "","",k_fit_min,k_fit_max);

    //k_fit.DrawCopy("Same");


    TF1 gauss_bckg("gauss_bckg", "[0]*TMath::Gaus(x,[1],[2])", rho_fit_min, 1400);
    gauss_bckg.SetParameters(7.81456e+05, -4.18595e+03, 1.37289e+03);
    rho_masses->ProjectionX()->Fit(&gauss_bckg, "","",900,1400);

    //gauss_bckg.SetLineColor(1);
    //gauss_bckg.DrawCopy("Same");


    TString rho_fit_string = create_interval(400., rho_fit_min)+"*("+get_par(linear_bckg, 0)+"+"+get_par(linear_bckg, 1)+"*x)+"
    +get_par(k_fit, 0)+"*TMath::Gaus(x,"+get_par(k_fit, 1)+","+get_par(k_fit, 2)+")+"
    +create_interval(rho_fit_min, 1400.)+"*("+get_par(gauss_bckg, 0)+"*TMath::Gaus(x,"+get_par(gauss_bckg, 1)+","+get_par(gauss_bckg, 2)
    +"))+"+create_interval(200., 1400.)+"*[0]*TMath::Gaus(x,[1],[2])";
    cout << rho_fit_string << endl;

    TF1 rho_fit("rho_fit", rho_fit_string, 200, 1400);
    rho_fit.SetParameters(1.93039e+03, rho_mass, 1.17468e+02);

    rho_masses->ProjectionX()->Fit(&rho_fit, "","",rho_fit_min,rho_fit_max);

    rho_fit.DrawCopy("Same");
    


    TString true_rho_string = get_par(rho_fit, 0)+"*TMath::Gaus(x,"+get_par(rho_fit, 1)+","+get_par(rho_fit, 2)+")";
    TF1 true_rho("true_rho", true_rho_string, 200, 1400);
    true_rho.SetLineColor(4);
    true_rho.DrawCopy("Same");
    
    


    TString true_k_string = get_par(k_fit, 0)+"*TMath::Gaus(x,"+get_par(k_fit, 1)+","+get_par(k_fit, 2)+")";
    TF1 true_k("true_k", true_k_string, 200, 1400);
    true_k.SetLineColor(3);
    true_k.DrawCopy("Same");


    TLatex text3(750,800,"\\rho");
    text3.DrawClone();

    TLatex text4(500,950,"K^{0}_{s}");
    text4.DrawClone();



    ///////////////////////////////////////////////////////


    projections->cd(2);
    rho_masses->ProjectionY()->SetTitle("Mass of second potential rho particle");
    rho_masses->ProjectionY()->GetYaxis()->SetTitle("events");
    rho_masses->ProjectionY()->Draw();
    //TLine line4 = TLine(rho_mass, 0, rho_mass, 1.05*rho_masses->ProjectionY()->GetMaximum());
    //line4.DrawClone();


    const float rho_fit_min_2 = 606.175;
    const float rho_fit_max_2 = 850.;
    const float k_fit_min_2 = 470.;
    const float k_fit_max_2 = 530.;

    TString linear_bckg_string_2 = "("+create_interval(k_fit_max_2, rho_fit_min_2)+"+"+create_interval(0., k_fit_min_2)+")*([0]+[1]*x)";
    cout << linear_bckg_string_2 << endl;
    TF1 linear_bckg_2("linear_bckg_2", linear_bckg_string_2, 200, 1400);
    linear_bckg_2.SetParameters(1, 1);

    rho_masses->ProjectionY()->Fit(&linear_bckg_2, "","",400,rho_fit_min_2);


    TString k_fit_string_2 = get_par(linear_bckg_2, 0)+"+"+get_par(linear_bckg_2, 1)+"*x+[0]*TMath::Gaus(x,[1],[2])";
    cout << k_fit_string_2 << endl;
    TF1 k_fit_2("k_fit", k_fit_string_2, 200, 1400);
    k_fit_2.SetParameters(8.19813e+02, 4.97588e+02, 1.19288e+01);

    rho_masses->ProjectionY()->Fit(&k_fit_2, "","",425,575);

    //k_fit.DrawCopy("Same");


    TF1 gauss_bckg_2("gauss_bckg_2", "[0]*TMath::Gaus(x,[1],[2])", rho_fit_min_2, 1400);
    gauss_bckg_2.SetParameters(1.22029e+07, -5.53468e+03, 1.45704e+03);
    rho_masses->ProjectionY()->Fit(&gauss_bckg_2, "","",900,1400);

    //gauss_bckg_2.SetLineColor(1);
    //gauss_bckg_2.DrawCopy("Same");


    TString rho_fit_string_2 = create_interval(400., rho_fit_min_2)+"*("+get_par(linear_bckg_2, 0)+"+"+get_par(linear_bckg_2, 1)+"*x)+"
    +get_par(k_fit_2, 0)+"*TMath::Gaus(x,"+get_par(k_fit_2, 1)+","+get_par(k_fit_2, 2)+")+"
    +create_interval(rho_fit_min_2, 1400.)+"*("+get_par(gauss_bckg_2, 0)+"*TMath::Gaus(x,"+get_par(gauss_bckg_2, 1)+","+get_par(gauss_bckg_2, 2)+"))+"
    +create_interval(rho_fit_min_2, 1400.)+"*[0]*TMath::Gaus(x,[1],[2])";
    cout << rho_fit_string_2 << endl;

    TF1 rho_fit_2("rho_fit_2", rho_fit_string_2, 200, 1400);
    rho_fit_2.SetParameters(1.93039e+03, rho_mass, 1.17468e+02);

    rho_masses->ProjectionY()->Fit(&rho_fit_2, "","",rho_fit_min_2,rho_fit_max_2);

    rho_fit_2.DrawCopy("Same");



    TString true_rho_string_2 = get_par(rho_fit_2, 0)+"*TMath::Gaus(x,"+get_par(rho_fit_2, 1)+","+get_par(rho_fit_2, 2)+")";
    TF1 true_rho_2("true_rho_2", true_rho_string_2, 200, 1400);
    true_rho_2.SetLineColor(4);
    true_rho_2.DrawCopy("Same");


    TString true_k_string_2 = get_par(k_fit_2, 0)+"*TMath::Gaus(x,"+get_par(k_fit_2, 1)+","+get_par(k_fit_2, 2)+")";
    TF1 true_k_2("true_k_2", true_k_string_2, 200, 1400);
    true_k_2.SetLineColor(3);
    true_k_2.DrawCopy("Same");


    TLatex text5(750,800,"\\rho");
    text5.DrawClone();

    TLatex text6(500,950,"K^{0}_{s}");
    text6.DrawClone();


    projections->Print("projections.pdf","RECREATE");
    
}