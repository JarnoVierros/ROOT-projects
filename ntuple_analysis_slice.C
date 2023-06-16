

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
        delta  = enforce_interval(projection->GetBinCenter(i), 650, 850)*(projection->GetBinContent(i)-func(projection->GetBinCenter(i),par));
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    cout << "chisq: " << chisq << endl << endl;
}

void ntuple_analysis_slice() {

    const string filename = "TOTEM43.root"; //"TOTEM43.root", kpkm.roo, 110000.root
    const float muon_mass = 139.57039;
    const float rho_mass = 730; //770

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

    const int rho_x_min = 550;
    const int rho_x_max = 1000;
    const int rho_y_min = 550;
    const int rho_y_max = 1000;

    TLine rho_line1 = TLine(rho_x_min, rho_y_min, rho_x_max, rho_y_min);
    rho_line1.DrawClone();

    TLine rho_line2 = TLine(rho_x_max, rho_y_min, rho_x_max, rho_y_max);
    rho_line2.DrawClone();

    TLine rho_line3 = TLine(rho_x_min, rho_y_max, rho_x_max, rho_y_max);
    rho_line3.DrawClone();

    TLine rho_line4 = TLine(rho_x_min, rho_y_min, rho_x_min, rho_y_max);
    rho_line4.DrawClone();


    const int k1_x_min = 200;
    const int k1_x_max = 800;
    const int k1_y_min = 450;
    const int k1_y_max = 550;

    TLine k1_line1 = TLine(k1_x_min, k1_y_min, k1_x_max, k1_y_min);
    k1_line1.DrawClone();

    TLine k1_line2 = TLine(k1_x_max, k1_y_min, k1_x_max, k1_y_max);
    k1_line2.DrawClone();

    TLine k1_line3 = TLine(k1_x_min, k1_y_max, k1_x_max, k1_y_max);
    k1_line3.DrawClone();

    TLine k1_line4 = TLine(k1_x_min, k1_y_min, k1_x_min, k1_y_max);
    k1_line4.DrawClone();

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


    auto projections = new TCanvas("Canvas2","Canvas2");

    
    TCutG *rho_cut = new TCutG("rho_cut",4);
    rho_cut->SetPoint(0,rho_x_min,rho_y_min);
    rho_cut->SetPoint(1,rho_x_max,rho_y_min);
    rho_cut->SetPoint(2,rho_x_max,rho_y_max);
    rho_cut->SetPoint(3,rho_x_min,rho_y_max);
    rho_cut->SetPoint(4,rho_x_min,rho_y_min);

    TCutG *k1_cut = new TCutG("k1_cut",4);
    k1_cut->SetPoint(0,k1_x_min,k1_y_min);
    k1_cut->SetPoint(1,k1_x_max,k1_y_min);
    k1_cut->SetPoint(2,k1_x_max,k1_y_max);
    k1_cut->SetPoint(3,k1_x_min,k1_y_max);
    k1_cut->SetPoint(4,k1_x_min,k1_y_min);

    TCutG *k2_cut = new TCutG("k2_cut",4);
    k2_cut->SetPoint(0,k2_x_min,k2_y_min);
    k2_cut->SetPoint(1,k2_x_max,k2_y_min);
    k2_cut->SetPoint(2,k2_x_max,k2_y_max);
    k2_cut->SetPoint(3,k2_x_min,k2_y_max);
    k2_cut->SetPoint(4,k2_x_min,k2_y_min);


    //rho_masses->ProjectionX("X_projection",1,200,"[-k1_cut, -k2_cut]")->Draw();
    projection = rho_masses->ProjectionX("X_projection",1,200,"[-k1_cut, -k2_cut]");
    projection->Draw();

    TMinuit *gMinuit = new TMinuit(3);
    gMinuit->SetFCN(fcn);
    
    //minuit->SetFitObject(X_projection)

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
    
}