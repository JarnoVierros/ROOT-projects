
const float peak = 780;//722.311

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

        float get_squared_total_dxy() {
            float total_squared_dxy = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_squared_dxy += pow(current_particle->dxy, 2);
            }
            return total_squared_dxy;
        }

        float get_squared_total_dz() {
            float total_squared_dz = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_squared_dz += pow(current_particle->dz, 2);
            }
            return total_squared_dz;
        }
};

Double_t func_gaussian(float x,Double_t *par){
    Double_t value = par[0]*TMath::Gaus(x,par[1],par[2]);
    return value;
}

Double_t func_BreitWigner(float x,Double_t *par){
    Double_t value = par[0]*TMath::BreitWigner(x,par[1],par[2]);
    return value;
}

Double_t func_Landau(float x,Double_t *par){
    Double_t value = par[0]*TMath::Landau(x,par[1],par[2]);
    return value;
}


static TH1D* projection;
static TH1D* raw_projection;

void fcn_gaussian(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(projection->GetBinCenter(i), peak-100, peak+100)*(projection->GetBinContent(i)-func_gaussian(projection->GetBinCenter(i),par))/projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void raw_fcn_gaussian(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (raw_projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(raw_projection->GetBinCenter(i), peak-100, peak+100)*(raw_projection->GetBinContent(i)-func_gaussian(raw_projection->GetBinCenter(i),par))/raw_projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void fcn_BreitWigner(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(projection->GetBinCenter(i), peak-100, peak+100)*(projection->GetBinContent(i)-func_BreitWigner(projection->GetBinCenter(i),par))/projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void fcn_Landau(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(projection->GetBinCenter(i), peak-100, 1400)*(projection->GetBinContent(i)-func_Landau(projection->GetBinCenter(i),par))/(projection->GetBinError(i));
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void ntuple_analysis_momentum_conservation_v2() {

    const string filename = "./ntuples/rho.root"; //"TOTEM43.root", kpkm.roo, 110000.root, rho.root, MinBias.root
    const bool monte_carlo = true;
    const float pion_mass = 139.57039;
    const float rho_mass = 730; //770

    const float allowed_px_difference = 150;
    const float allowed_py_difference = 150;
    const float allowed_dxy_variance = 0.1;
    const float allowed_dz_variance = 0.6;
    const float allowed_rho_2_mass_difference = 150;

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

    auto total_dxy_squared_vs_rho_m1 = new TH2F("total_dxy_squared_vs_rho_m1", "total dxy^2 and rho1 mass;total dxy^2;rho1 mass/MeV",200,0,5,200,200,1400);
    auto total_dz_squared_vs_rho_m1 = new TH2F("total_dz_squared_vs_rho_m1", "total dz^2 and rho1 mass;total dz^2;rho1 mass/MeV",200,0,5,200,200,1400);

    auto total_dxy_squared_vs_rho_m2 = new TH2F("total_dxy_squared_vs_rho_m2", "total dxy^2 and rho2 mass;total dxy^2;rho2 mass/MeV",200,0,5,200,200,1400);
    auto total_dz_squared_vs_rho_m2 = new TH2F("total_dz_squared_vs_rho_m2", "total dz^2 and rho2 mass;total dz^2;rho2 mass/MeV",200,0,5,200,200,1400);

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

            total_dxy_squared_vs_rho_m1->Fill(current_event.get_squared_total_dxy(), raw_masses[0]);
            total_dz_squared_vs_rho_m1->Fill(current_event.get_squared_total_dz(), raw_masses[0]);

            total_dxy_squared_vs_rho_m2->Fill(current_event.get_squared_total_dxy(), raw_masses[1]);
            total_dz_squared_vs_rho_m2->Fill(current_event.get_squared_total_dz(), raw_masses[1]);
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
    raw_projection = rho_masses_raw->ProjectionX("X_projection_raw");
    raw_projection->Draw();

    float results[4];

    TMinuit *gMinuit0 = new TMinuit(3);
    gMinuit0->SetFCN(raw_fcn_gaussian);

    Double_t arglist0[10];
    Int_t ierflg0 = 0;

    arglist0[0] = 1;

    gMinuit0->mnexcm("SET ERR", arglist0 ,1,ierflg0);

    static Double_t vstart0[4] = {1.92522e+03, 7.14645e+02, 1.24676e+02};
    static Double_t step0[4] = {1, 1, 1};
    gMinuit0->mnparm(0, "a0", vstart0[0], step0[0], 0,0,ierflg0);
    gMinuit0->mnparm(1, "a2", vstart0[1], step0[1], 0,0,ierflg0);
    gMinuit0->mnparm(2, "a3", vstart0[2], step0[2], 0,0,ierflg0);

    arglist0[0] = 500;
    arglist0[1] = 1.;
    gMinuit0->mnexcm("MIGRAD", arglist0 ,2,ierflg0);

    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit0->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    Double_t val1,err1,val2,err2,val3,err3;
    gMinuit0->GetParameter(0, val1, err1);
    gMinuit0->GetParameter(1, val2, err2);
    results[0] = val2;
    gMinuit0->GetParameter(2, val3, err3);

    TF1 gaussia_fit_raw("gaussia_fit_raw", "[0]*TMath::Gaus(x,[1],[2])", 200, 1400);
    gaussia_fit_raw.SetParameters(val1, val2, val3);

    gaussia_fit_raw.DrawCopy("Same");


    auto main = new TCanvas("Canvas4","Canvas4");
    rho_masses->Draw("Colz");

 

    auto projections = new TCanvas("Canvas5","Canvas5");
    

    //projection = rho_masses->ProjectionX("X_projection");
    projection = rho_masses->ProjectionX("X_projection", (rho_mass-allowed_rho_2_mass_difference-200)/6, (rho_mass+allowed_rho_2_mass_difference-200)/6);
    projection->Draw();


/////////////////////////////////////

    TMinuit *gMinuit1 = new TMinuit(3);
    gMinuit1->SetFCN(fcn_gaussian);

    Double_t arglist1[10];
    Int_t ierflg1 = 0;

    arglist1[0] = 1;

    gMinuit1->mnexcm("SET ERR", arglist1 ,1,ierflg1);

    static Double_t vstart1[4] = {6.16160e+02, 7.22303e+02, 1.17525e+02};
    static Double_t step1[4] = {1, 1, 1};
    gMinuit1->mnparm(0, "a1", vstart1[0], step1[0], 0,0,ierflg1);
    gMinuit1->mnparm(1, "a2", vstart1[1], step1[1], 0,0,ierflg1);
    gMinuit1->mnparm(2, "a3", vstart1[2], step1[2], 0,0,ierflg1);

    arglist1[0] = 500;
    arglist1[1] = 1.;
    gMinuit1->mnexcm("MIGRAD", arglist1 ,2,ierflg1);

    gMinuit1->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    
    gMinuit1->GetParameter(0, val1, err1);
    gMinuit1->GetParameter(1, val2, err2);
    results[1] = val2;
    gMinuit1->GetParameter(2, val3, err3);

    TF1 gaussia_fit("gaussia_fit", "[0]*TMath::Gaus(x,[1],[2])", 200, 1400);
    gaussia_fit.SetParameters(val1, val2, val3);

    gaussia_fit.DrawCopy("Same");

    
////////////////////////////////////////////
    
    TMinuit *gMinuit2 = new TMinuit(3);
    gMinuit2->SetFCN(fcn_BreitWigner);

    Double_t arglist2[10];
    Int_t ierflg2 = 0;

    arglist2[0] = 1;

    gMinuit2 ->mnexcm("SET ERR", arglist2 ,1,ierflg2);

    static Double_t vstart2[4] = {3.26161e+05, 7.11412e+02, 3.42154e+02};
    static Double_t step2[4] = {1, 1, 1};
    gMinuit2 ->mnparm(0, "a1", vstart2[0], step2[0], 0,0,ierflg2);
    gMinuit2 ->mnparm(1, "a2", vstart2[1], step2[1], 0,0,ierflg2);
    gMinuit2 ->mnparm(2, "a3", vstart2[2], step2[2], 0,0,ierflg2);

    arglist2[0] = 500;
    arglist2[1] = 1.;
    gMinuit2 ->mnexcm("MIGRAD", arglist2 ,2,ierflg2);

    //Double_t amin,edm,errdef;
    //Int_t nvpar,nparx,icstat;
    gMinuit2 ->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    //Double_t val1,err1,val2,err2,val3,err3;
    gMinuit2 ->GetParameter(0, val1, err1);
    gMinuit2 ->GetParameter(1, val2, err2);
    results[2] = val2;
    gMinuit2 ->GetParameter(2, val3, err3);

    TF1 BreitWigner_fit("BreitWigner_fit", "[0]*TMath::BreitWigner(x,[1],[2])", 200, 1400);
    BreitWigner_fit.SetParameters(val1, val2, val3);
    BreitWigner_fit.SetLineColor(3);
    BreitWigner_fit.DrawCopy("Same");


////////////////////////////////////////////

    TMinuit *gMinuit3 = new TMinuit(3);
    gMinuit3->SetFCN(fcn_Landau);

    Double_t arglist3[10];
    Int_t ierflg3 = 0;

    arglist3[0] = 1;

    gMinuit3 ->mnexcm("SET ERR", arglist3 ,1,ierflg3);

    static Double_t vstart3[4] = {3.28034e+03, 7.12024e+02, 7.82408e+01};
    static Double_t step3[4] = {1, 0.1, 1};
    gMinuit3 ->mnparm(0, "a1", vstart3[0], step3[0], 0,0,ierflg3);
    gMinuit3 ->mnparm(1, "a2", vstart3[1], step3[1], 0,0,ierflg3);
    gMinuit3 ->mnparm(2, "a3", vstart3[2], step3[2], 0,0,ierflg3);

    arglist3[0] = 500;
    arglist3[1] = 1.;
    gMinuit3 ->mnexcm("MIGRAD", arglist3 ,2,ierflg3);

    //Double_t amin,edm,errdef;
    //Int_t nvpar,nparx,icstat;
    gMinuit3 ->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    gMinuit3 ->GetParameter(0, val1, err1);
    gMinuit3 ->GetParameter(1, val2, err2);
    results[3] = val2;
    gMinuit3 ->GetParameter(2, val3, err3);

    TF1 Landau_fit("Landau_fit", "[0]*TMath::Landau(x,[1],[2])", 200, 1400);
    Landau_fit.SetParameters(val1, val2, val3);
    Landau_fit.SetLineColor(4);
    Landau_fit.DrawCopy("Same");



    

    auto dxy_comparison_1 = new TCanvas("Canvas6","Canvas6");
    total_dxy_squared_vs_rho_m1->Draw("Colz");

    auto dz_comparison_1 = new TCanvas("Canvas7","Canvas7");
    total_dz_squared_vs_rho_m1->Draw("Colz");

    auto dxy_comparison_2 = new TCanvas("Canvas8","Canvas8");
    total_dxy_squared_vs_rho_m2->Draw("Colz");

    auto dz_comparison_2 = new TCanvas("Canvas9","Canvas9");
    total_dz_squared_vs_rho_m2->Draw("Colz");

    cout << "RESULTS: " << endl;
    cout << results[0] << endl;
    cout << results[1] << endl;
    cout << results[2] << endl;
    cout << results[3] << endl;

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

/*
new TOTEM43
710.461
720.492
721.018
713.57
*/