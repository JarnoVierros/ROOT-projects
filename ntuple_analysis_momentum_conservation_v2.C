

string create_interval(float start, float stop) {
    return "((("+to_string(start)+"<x) ? 1 : 0) - (("+to_string(stop)+"<x) ? 1 : 0))";
}

int enforce_interval(float x, float start, float stop) {
    return ((start<x ? 1 : 0) - (stop<x ? 1 : 0));
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

        float dxy;
        float dz;

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
        const float prot_momentum = 6.5e+6;

        Particle** particles;
        int particle_count;

        float ThxR;
        float ThyR;
        float ThxL;
        float ThyL;

        float PR_p[3] = {0,0,prot_momentum};
        float PL_p[3] = {0,0,prot_momentum};


        float ref_p[3] = {0,0,0}; //momentums of the refractive system: 0=p_x, 1=p_y, 2=p_z

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

        void calculate_proton_momentums() {
            PR_p[0] = ThxR*prot_momentum;
            PL_p[0] = ThxL*prot_momentum;
            PR_p[1] = ThyR*prot_momentum;
            PL_p[1] = ThyL*prot_momentum;
        }

        void calculate_momentum_of_refractive_system() {
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                ref_p[0] += current_particle->p_x;
                ref_p[1] += current_particle->p_y;
                ref_p[2] += current_particle->p_z;
            }
        }

        float get_total_dxy() {
            float total_dxy = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_dxy += current_particle->dxy;
            }
            return total_dxy;
        }

        float get_total_dz() {
            float total_dz = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_dz += current_particle->dz;
            }
            return total_dz;
        }

        float get_dxy_variance() {
            float average_dxy = get_total_dxy()/particle_count;
            float variance = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                variance += pow(average_dxy-current_particle->dxy, 2);
            }
            return variance;
        }

        float get_dz_variance() {
            float average_dz = get_total_dz()/particle_count;
            float variance = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                variance += pow(average_dz-current_particle->dz, 2);
            }
            return variance;
        }
};

Double_t func(float x,Double_t *par){
    Double_t value = par[0]*TMath::Gaus(x,par[1],par[2]);
    /*
    cout << "x: " << x << endl;
    cout << "par[0]: " << par[0] << endl;
    cout << "par[1]: " << par[1] << endl;
    cout << "par[2]: " << par[2] << endl;
    cout << "VALUE: " << value << endl << endl;
    */
    return value;
}


static TH1D* projection;
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(projection->GetBinCenter(i), 650, 850)*(projection->GetBinContent(i)-func(projection->GetBinCenter(i),par))/projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    cout << "chisq: " << chisq << endl << endl;
}

void ntuple_analysis_momentum_conservation_v2() {

    const string filename = "TOTEM43.root"; //"TOTEM43.root", kpkm.roo, 110000.root, rho.root, MinBias.root
    const bool monte_carlo = false;
    const float pion_mass = 139.57039;
    const float rho_mass = 730; //770

    const float allowed_px_difference = 150;
    const float allowed_py_difference = 150;
    const float allowed_dxy_variance = 0.1;
    const float allowed_dz_variance = 0.1;
    const float allowed_rho_2_mass_difference = 100;

    TFile *file = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)file->Get("tree");
    TTreeReader Reader(tree);
    //TTreeReader Reader("tree", file);

    auto rho_masses = new TH2F("rho_masses", "Masses of potential rho particles;m1/MeV;m2/MeV",200,200,1400,200,200,1400);
    auto rho_masses_raw = new TH2F("rho_masses_raw", "Masses of raw potential rho particles;m1/MeV;m2/MeV",200,200,1400,200,200,1400);

    auto prot_px_vs_ref_px = new TH2F("prot_px_vs_ref_px", "x-momentum of protons and refractive system;prot_px/MeV;ref_px/MeV",200,-2000,2000,200,-2000,2000);
    auto prot_py_vs_ref_py = new TH2F("prot_py_vs_ref_py", "y-momentum of protons and refractive system;prot_py/MeV;ref_py/MeV",200,-2000,2000,200,-2000,2000);
    

    //auto total_dxy_vs_rho_m = new TH2F("total_dxy_vs_rho_m", "total dxy and rho mass;total dxy;rho mass/MeV",200,-5,5,200,200,1400);
    //auto total_dz_vs_rho_m = new TH2F("total_dz_vs_rho_m", "total dz and rho mass;total dz;rho mass/MeV",200,-5,5,200,200,1400);

    auto dxy_variance_vs_rho_m1 = new TH2F("dxy_variance_vs_rho_m1", "dxy variance and rho1 mass;dxy variance;rho1 mass/MeV",200,0,5,200,200,1400);
    auto dz_variance_vs_rho_m1 = new TH2F("dz_variance_vs_rho_m1", "dz variance and rho1 mass;dz variance;rho1 mass/MeV",200,0,5,200,200,1400);

    auto dxy_variance_vs_rho_m2 = new TH2F("dxy_variance_vs_rho_m2", "dxy variance and rho2 mass;dxy variance;rho2 mass/MeV",200,0,5,200,200,1400);
    auto dz_variance_vs_rho_m2 = new TH2F("dz_variance_vs_rho_m2", "dz variance and rho2 mass;dz variance;rho2 mass/MeV",200,0,5,200,200,1400);

    TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
    TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");
    TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
    TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
    TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");

    TTreeReaderArray<Float_t> trk_dxy(Reader, "trk_dxy");
    TTreeReaderArray<Float_t> trk_dz(Reader, "trk_dz");

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
            particle->m = pion_mass;
            particle->calculate_3d_momentum();
            particle->calculate_energy();
            particles[i] = particle;

            particle->dxy = trk_dxy[i];
            particle->dz = trk_dz[i];
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

        if (monte_carlo) {
            current_event.ThxR = -(*ThxR);
            current_event.ThyR = -(*ThyR);
            current_event.ThxL = (*ThxL);
            current_event.ThyL = (*ThyL);
        } else {
            current_event.ThxR = -(*ThxR);
            current_event.ThyR = (*ThyR);
            current_event.ThxL = -(*ThxL);
            current_event.ThyL = (*ThyL);
        }
        

        float total_charge = current_event.calculate_total_charge();
        if (total_charge != 0) {
            //cout << "INVALID" << endl << endl;
            continue;
        }

        current_event.calculate_momentum_of_refractive_system();
        current_event.calculate_proton_momentums();

        /*
        for (int i=0; i<3; ++i) {
            cout << current_event.ref_p[i] << endl;
        }

        cout << endl;

        cout << current_event.ThxR << endl;
        cout << current_event.ThyR << endl;
        for (int i=0; i<3; ++i) {
            cout << current_event.PR_p[i] << endl;
        }

        cout << endl;

        cout << current_event.ThxL << endl;
        cout << current_event.ThyL << endl;
        for (int i=0; i<3; ++i) {
            cout << current_event.PL_p[i] << endl;
        }
        */

        Particle* rhos[2][2];
        current_event.reconstruct_2_rhos(rhos);

        prot_px_vs_ref_px->Fill(current_event.PR_p[0]+current_event.PL_p[0], current_event.ref_p[0]);
        prot_py_vs_ref_py->Fill(current_event.PR_p[1]+current_event.PL_p[1], current_event.ref_p[1]);
        

        float raw_masses[2];
        for (int i=0; i<2; ++i) {
            for (int j=0; j<2; ++j) {
                Particle* rho = rhos[i][j];
                raw_masses[j] = rho->m;
                //cout << "m: " << rho->m << endl;
            }
            rho_masses_raw->Fill(raw_masses[0], raw_masses[1]);

            dxy_variance_vs_rho_m1->Fill(current_event.get_dxy_variance(), raw_masses[0]);
            dz_variance_vs_rho_m1->Fill(current_event.get_dz_variance(), raw_masses[0]);

            dxy_variance_vs_rho_m2->Fill(current_event.get_dxy_variance(), raw_masses[1]);
            dz_variance_vs_rho_m2->Fill(current_event.get_dz_variance(), raw_masses[1]);
        }
        
        if (allowed_px_difference < abs(abs(current_event.PR_p[0]+current_event.PL_p[0]) - abs(current_event.ref_p[0]))) {
            continue;
        }

        if (allowed_py_difference < abs(abs(current_event.PR_p[1]+current_event.PL_p[1]) - abs(current_event.ref_p[1]))) {
            continue;
        }
        
        if (current_event.get_dxy_variance()>allowed_dxy_variance) {
            continue;
        }

        if (current_event.get_dz_variance()>allowed_dz_variance) {
            continue;
        }

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
        //cout << endl;

    } 


    auto px_comparison = new TCanvas("Canvas0","Canvas0");
    prot_px_vs_ref_px->Draw("Colz");

    TLine upper_limit_x = TLine(-1500, 1500+allowed_px_difference, 1500, -1500+allowed_px_difference);
    upper_limit_x.DrawClone();

    TLine lower_limit_x = TLine(-1500, 1500-allowed_px_difference, 1500, -1500-allowed_px_difference);
    lower_limit_x.DrawClone();


    auto py_comparison = new TCanvas("Canvas1","Canvas1");
    prot_py_vs_ref_py->Draw("Colz");

    TLine upper_limit_y = TLine(-1500, 1500+allowed_py_difference, 1500, -1500+allowed_py_difference);
    upper_limit_y.DrawClone();

    TLine lower_limit_y = TLine(-1500, 1500-allowed_py_difference, 1500, -1500-allowed_py_difference);
    lower_limit_y.DrawClone();


    auto raw = new TCanvas("Canvas2","Canvas2");
    rho_masses_raw->Draw("Colz");
    //rho_masses_raw->ProjectionX("X_projection_raw")->Draw();

    auto raw_projections = new TCanvas("Canvas3","Canvas3");
    rho_masses_raw->ProjectionX("X_projection_raw")->Draw();

    auto main = new TCanvas("Canvas4","Canvas4");
    rho_masses->Draw("Colz");

    /*
    const int k2_x_min = 450;
    const int k2_x_max = 550;
    const int k2_y_min = 200;
    const int k2_y_max = 800;

    TLine k2_line1 = TLine(k2_x_min, k2_y_min, k2_x_max, k2_y_min);
    k2_line1.DrawClone();

    TLine k2_line2 = TLine(k2_x_max, k2_y_min, k2_x_max, k2_y_max);
    k2_line2.DrawClone();

    TLine k2_line3 = TLine(k2_x_min, k2_y_max, k2_x_max, k2_y_max);
    k2_line3.DrawClone();

    TLine k2_line4 = TLine(k2_x_min, k2_y_min, k2_x_min, k2_y_max);
    k2_line4.DrawClone();


    TCutG *k2_cut = new TCutG("k2_cut",4);
    k2_cut->SetPoint(0,k2_x_min,k2_y_min);
    k2_cut->SetPoint(1,k2_x_max,k2_y_min);
    k2_cut->SetPoint(2,k2_x_max,k2_y_max);
    k2_cut->SetPoint(3,k2_x_min,k2_y_max);
    k2_cut->SetPoint(4,k2_x_min,k2_y_min);


    TLine ycut_line1 = TLine(200, 600, 1400, 600);
    ycut_line1.DrawClone();

    TLine ycut_line2 = TLine(200, 800, 1400, 800);
    ycut_line2.DrawClone();
    */

    auto projections = new TCanvas("Canvas5","Canvas5");
    

    //projection = rho_masses->ProjectionX("X_projection");
    projection = rho_masses->ProjectionX("X_projection", (rho_mass-allowed_rho_2_mass_difference-200)/6, (rho_mass+allowed_rho_2_mass_difference-200)/6);
    projection->Draw();


    

    auto dxy_comparison_1 = new TCanvas("Canvas6","Canvas6");
    dxy_variance_vs_rho_m1->Draw("Colz");

    auto dz_comparison_1 = new TCanvas("Canvas7","Canvas7");
    dz_variance_vs_rho_m1->Draw("Colz");

    auto dxy_comparison_2 = new TCanvas("Canvas8","Canvas8");
    dxy_variance_vs_rho_m2->Draw("Colz");

    auto dz_comparison_2 = new TCanvas("Canvas9","Canvas9");
    dz_variance_vs_rho_m2->Draw("Colz");

    /*
    TMinuit *gMinuit = new TMinuit(3);
    gMinuit->SetFCN(fcn);

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;

    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    static Double_t vstart[4] = {1927.99, 730, 118.995};
    static Double_t step[4] = {1, 1, 1};
    gMinuit->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);

    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    cout << "amin: " << amin << endl;
    cout << "edm: " << edm << endl;
    cout << "errdef: " << errdef << endl;
    cout << "nvpar: " << nvpar << endl;
    cout << "nparx: " << nparx << endl;
    cout << "icstat: " << icstat << endl;
    Double_t val1,err1,val2,err2,val3,err3;
    gMinuit->GetParameter(0, val1, err1);
    gMinuit->GetParameter(1, val2, err2);
    gMinuit->GetParameter(2, val3, err3);
    cout << "parameter 1 val/err: " << val1 << " / " << err1 << endl;
    cout << "parameter 2 val/err: " << val2 << " / " << err2 << endl;
    cout << "parameter 3 val/err: " << val3 << " / " << err3 << endl;

    TF1 gaussia_fit("gaussia_fit", "[0]*TMath::Gaus(x,[1],[2])", 200, 1400);
    gaussia_fit.SetParameters(val1, val2, val3);

    gaussia_fit.DrawCopy("Same");
    */   
}