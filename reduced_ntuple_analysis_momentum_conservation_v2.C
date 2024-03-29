
const float global_peak = 780;//722.311

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
        //unknown=0, pion=1, rho=2, glueball=3
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

        void reconstruct_2_from_4(Particle* origins[2][2], int type) {
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

            Particle* origin1_c1 = reconstruct_1_from_2(positives[0], negatives[0], type);
            origins[0][0] = origin1_c1;
            Particle* origin2_c1 = reconstruct_1_from_2(positives[1], negatives[1], type);
            origins[0][1] = origin2_c1;
            
            Particle* origin1_c2 = reconstruct_1_from_2(positives[0], negatives[1], type);
            origins[1][0] = origin1_c2;
            Particle* origin2_c2 = reconstruct_1_from_2(positives[1], negatives[0], type);
            origins[1][1] = origin2_c2;
        }
    
        Particle* reconstruct_1_from_2(Particle* particle1, Particle* particle2, int type) {
            Particle* origin = new Particle(type);
            origin->E = particle1->E + particle2->E;
            origin->p_x = particle1->p_x + particle2->p_x;
            origin->p_y = particle1->p_y + particle2->p_y;
            origin->p_z = particle1->p_z + particle2->p_z;
            origin->calculate_total_momentum();
            origin->calculate_mass();
            return origin;
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
    Double_t value = par[3] + par[0]*TMath::Landau(x,par[1],par[2]);
    return value;
}


static TH1D* projection;
static TH1D* raw_projection;

static float fcn_gaussian_min;
static float fcn_gaussian_max;
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
        //cout << "HERE1" << endl;
        //cout << fcn_gaussian_min << endl;
        //cout << fcn_gaussian_max << endl;
        delta  = enforce_interval(projection->GetBinCenter(i), fcn_gaussian_min, fcn_gaussian_max)*(projection->GetBinContent(i)-func_gaussian(projection->GetBinCenter(i),par))/projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

static float raw_fcn_gaussian_min;
static float raw_fcn_gaussian_max;
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
        //cout << "HERE2" << endl;
        //cout << raw_fcn_gaussian_min << endl;
        //cout << raw_fcn_gaussian_max << endl;
        delta  = enforce_interval(raw_projection->GetBinCenter(i), raw_fcn_gaussian_min, raw_fcn_gaussian_max)*(raw_projection->GetBinContent(i)-func_gaussian(raw_projection->GetBinCenter(i),par))/raw_projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

static float fcn_BreitWigner_min;
static float fcn_BreitWigner_max;
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
        delta  = enforce_interval(projection->GetBinCenter(i), fcn_BreitWigner_min, fcn_BreitWigner_max)*(projection->GetBinContent(i)-func_BreitWigner(projection->GetBinCenter(i),par))/projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

static float fcn_Landau_min;
static float fcn_Landau_max;
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
        delta  = enforce_interval(projection->GetBinCenter(i), fcn_Landau_min, fcn_Landau_max)*(projection->GetBinContent(i)-func_Landau(projection->GetBinCenter(i),par))/(projection->GetBinError(i));
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void reduced_ntuple_analysis_momentum_conservation_v2() {

    const string filename = "./ntuples/TOTEM.root"; //"TOTEM43.root", kpkm.roo, 110000.root, rho.root, MinBias.root
    const bool monte_carlo = true;

    const float pion_mass = 139.57039;
    const float rho_mass = 720; //770

    const float allowed_px_difference = 150;
    const float allowed_py_difference = 150;
    const float allowed_squared_total_dxy = 0.15;
    const float allowed_squared_total_dz = 0.2;
    const float allowed_rho_2_mass_difference = 150;
    const float allowed_rho_2_mass_difference_supercut = 31;

    TFile *file = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)file->Get("tree");
    TTreeReader Reader(tree);
    //TTreeReader Reader("tree", file);

    auto rho_masses = new TH2F("rho_masses", "Masses of potential rho particles with impact parameter cuts;m1/MeV;m2/MeV",200,200,1400,200,200,1400);
    auto rho_masses_raw = new TH2F("rho_masses_raw", "Masses of potential rho particles;m1/MeV;m2/MeV",200,200,1400,200,200,1400);

    auto raw_origin_mass = new TH1F("raw_origin_mass", "Mass of the four track origin without cuts;mass/MeV",200,500,3000);
    auto origin_mass = new TH1F("origin_mass", "Mass of the four track origin with cuts;mass/MeV",200,500,3000);

    auto raw_origin_m_vs_rho1_m = new TH2F("raw_origin_m_vs_rho1_m", "Mass of the four track origin compared with mass of first rho particle without cuts;origin mass/MeV;rho1 mass/MeV",200,500,3000,200,200,1400);
    auto origin_m_vs_rho1_m = new TH2F("origin_m_vs_rho1_m", "Mass of the four track origin compared with mass of first rho particle with cuts;origin mass/MeV;rho1 mass/MeV",200,500,3000,200,200,1400);

    auto origin_m_vs_rho1_m_supercut = new TH2F("origin_m_vs_rho1_m_supercut", "Mass of the four track origin compared with mass of first rho particle with supercuts;origin mass/MeV;rho1 mass/MeV",200,1000,3000,200,200,1400);

    origin_m_vs_rho1_m_supercut->Sumw2();

    TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
    TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");
    TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
    TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
    TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");

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
        }

        Event current_event(particles, 4);

        float total_charge = current_event.calculate_total_charge();
        if (total_charge != 0) {
            continue;
        }

        current_event.calculate_momentum_of_refractive_system();

        Particle* rhos[2][2];
        current_event.reconstruct_2_from_4(rhos, 2);

        Particle* glueball1 = current_event.reconstruct_1_from_2(rhos[0][0], rhos[0][1], 3);
        Particle* glueball2 = current_event.reconstruct_1_from_2(rhos[1][0], rhos[1][1], 3);

        raw_origin_mass->Fill(glueball1->m);
        raw_origin_mass->Fill(glueball2->m);

        raw_origin_m_vs_rho1_m->Fill(glueball1->m, rhos[0][0]->m);
        raw_origin_m_vs_rho1_m->Fill(glueball2->m, rhos[1][0]->m);
        

        float raw_masses[2];
        for (int i=0; i<2; ++i) {
            for (int j=0; j<2; ++j) {
                Particle* rho = rhos[i][j];
                raw_masses[j] = rho->m;
            }
            rho_masses_raw->Fill(raw_masses[0], raw_masses[1]);
        }

        float masses[2];
        for (int i=0; i<2; ++i) {
            for (int j=0; j<2; ++j) {
                Particle* rho = rhos[i][j];
                masses[j] = rho->m;
            }
            rho_masses->Fill(masses[0], masses[1]);
        }

        if (rho_mass - allowed_rho_2_mass_difference < rhos[0][1]->m && rhos[0][1]->m < rho_mass + allowed_rho_2_mass_difference ) {
            origin_mass->Fill(glueball1->m);
            origin_m_vs_rho1_m->Fill(glueball1->m, rhos[0][0]->m);
        }

        if (rho_mass - allowed_rho_2_mass_difference < rhos[1][1]->m && rhos[1][1]->m < rho_mass + allowed_rho_2_mass_difference ) {
            origin_mass->Fill(glueball2->m);
            origin_m_vs_rho1_m->Fill(glueball2->m, rhos[1][0]->m);
        }


        if (rho_mass - allowed_rho_2_mass_difference_supercut < rhos[0][1]->m && rhos[0][1]->m < rho_mass + allowed_rho_2_mass_difference_supercut) {
            origin_m_vs_rho1_m_supercut->Fill(glueball1->m, rhos[0][0]->m);
        }

        if (rho_mass - allowed_rho_2_mass_difference_supercut < rhos[1][1]->m && rhos[1][1]->m < rho_mass + allowed_rho_2_mass_difference_supercut) {
            origin_m_vs_rho1_m_supercut->Fill(glueball2->m, rhos[1][0]->m);
        }

        //cout << endl;
        //cout << endl;

    } 

    auto raw = new TCanvas("Canvas2","Canvas2");
    rho_masses_raw->Draw("Colz");

    auto raw_projections = new TCanvas("Canvas3","Canvas3");
    raw_projection = rho_masses_raw->ProjectionX("X_projection_raw");
    raw_projection->Draw();

    float results[4];

    int peak_bin = raw_projection->GetMaximumBin();
    float peak = raw_projection->GetBinCenter(peak_bin);

    raw_fcn_gaussian_min = peak - 100;
    raw_fcn_gaussian_max = peak + 100;

    TMinuit *gMinuit0 = new TMinuit(3);
    gMinuit0->SetFCN(raw_fcn_gaussian);

    Double_t arglist0[10];
    Int_t ierflg0 = 0;

    arglist0[0] = 1;

    gMinuit0->mnexcm("SET ERR", arglist0 ,1,ierflg0);

    static Double_t vstart0[4] = {1.92522e+03, peak, 1.24676e+02};
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

    Double_t val1,err1,val2,err2,val3,err3,val4,err4;
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

    peak_bin = projection->GetMaximumBin();
    peak = projection->GetBinCenter(peak_bin);

    float dynamic_rho_mass = peak;

    fcn_gaussian_min = peak - 100;
    fcn_gaussian_max = peak + 100;

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

    fcn_BreitWigner_min = peak - 100;;
    fcn_BreitWigner_max = peak + 100;;
    
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

    fcn_Landau_min = peak - 100;
    fcn_Landau_max = 1400;

    TMinuit *gMinuit3 = new TMinuit(4);
    gMinuit3->SetFCN(fcn_Landau);

    Double_t arglist3[10];
    Int_t ierflg3 = 0;

    arglist3[0] = 1;

    gMinuit3 ->mnexcm("SET ERR", arglist3 ,1,ierflg3);

    static Double_t vstart3[4] = {3.28034e+03, 7.12024e+02, 7.82408e+01, 0};
    static Double_t step3[4] = {1, 0.1, 1, 1};
    gMinuit3 ->mnparm(0, "a1", vstart3[0], step3[0], 0,0,ierflg3);
    gMinuit3 ->mnparm(1, "a2", vstart3[1], step3[1], 0,0,ierflg3);
    gMinuit3 ->mnparm(2, "a3", vstart3[2], step3[2], 0,0,ierflg3);
    gMinuit3 ->mnparm(3, "a4", vstart3[3], step3[3], 0,0,ierflg3);

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
    gMinuit3 ->GetParameter(3, val4, err4);

    TF1 Landau_fit("Landau_fit", "[3] + [0]*TMath::Landau(x,[1],[2])", 200, 1400);
    Landau_fit.SetParameters(val1, val2, val3, val4);
    Landau_fit.SetLineColor(4);
    Landau_fit.DrawCopy("Same");


    cout << "RESULTS: " << endl;
    cout << results[0] << endl;
    cout << results[1] << endl;
    cout << results[2] << endl;
    cout << results[3] << endl;


    auto raw_origin_mass_canvas = new TCanvas("Canvas10","Canvas10");
    raw_origin_mass->Draw("Colz");

    peak_bin = raw_origin_mass->GetMaximumBin();
    peak = raw_origin_mass->GetBinCenter(peak_bin);
    TF1 raw_origin_mass_fit("raw_origin_mass_fit", "[0]*TMath::Gaus(x,[1],[2])", 500, 3000);
    raw_origin_mass_fit.SetParameters(2000, peak, 100);
    //raw_origin_mass->Fit(&raw_origin_mass_fit, "","",peak-50,peak+50);


    auto origin_mass_canvas = new TCanvas("Canvas11","Canvas11");
    origin_mass->Draw("Colz");

    peak_bin = origin_mass->GetMaximumBin();
    peak = origin_mass->GetBinCenter(peak_bin);
    TF1 origin_mass_fit("origin_mass_fit", "[0]*TMath::Gaus(x,[1],[2])", 500, 3000);
    origin_mass_fit.SetParameters(2000, peak, 100);
    //origin_mass->Fit(&origin_mass_fit, "","",peak-50,peak+50);
    
    auto raw_origin_m_vs_rho1_m_canvas = new TCanvas("Canvas12","Canvas12");
    raw_origin_m_vs_rho1_m->Draw("Colz");

    auto origin_mass_canvas_canvas = new TCanvas("Canvas13","Canvas13");
    origin_m_vs_rho1_m->Draw("Colz");

    auto origin_mass_canvas_canvas_supercut = new TCanvas("Canvas14","Canvas14");
    origin_m_vs_rho1_m_supercut->Draw("Colz");

    TLine line1 = TLine(1000, dynamic_rho_mass-allowed_rho_2_mass_difference_supercut, 3000, dynamic_rho_mass-allowed_rho_2_mass_difference_supercut);
    line1.DrawClone();

    TLine line2 = TLine(1000, dynamic_rho_mass+allowed_rho_2_mass_difference_supercut, 3000, dynamic_rho_mass+allowed_rho_2_mass_difference_supercut);
    line2.DrawClone();

    auto origin_projection_supercut_canvas = new TCanvas("Canvas15","Canvas15");
    
    auto origin_projection_supercut = origin_m_vs_rho1_m_supercut->ProjectionX("origin_m_vs_rho1_m_projection_supercut", (dynamic_rho_mass-allowed_rho_2_mass_difference_supercut-200)/6, (dynamic_rho_mass+allowed_rho_2_mass_difference_supercut-200)/6);
    origin_projection_supercut->Draw();

}

/*
new TOTEM43
715.048
723.566
729.022
726.99
*/