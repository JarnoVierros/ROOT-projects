
#include "CMS_lumi.C"

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


        float diff_p[3] = {0,0,0}; //momentums of the diffractive system: 0=p_x, 1=p_y, 2=p_z

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

        void calculate_momentum_of_diffractive_system() {
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                diff_p[0] += current_particle->p_x;
                diff_p[1] += current_particle->p_y;
                diff_p[2] += current_particle->p_z;
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

        float get_greatest_dxy() {
            float maximum = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                float dxy = abs(current_particle->dxy);
                if (maximum < dxy) {
                    maximum = dxy;
                }
            }
            return maximum;
        }

        float get_greatest_dz() {
            float maximum = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                float dz = abs(current_particle->dz);
                if (maximum < dz) {
                    maximum = dz;
                }
            }
            return maximum;
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

static Double_t LanBreWig_parameters[4];
Double_t func_LanBreWig(float x,Double_t *par){
    Double_t value = par[0]*TMath::BreitWigner(x,par[1],par[2]) + LanBreWig_parameters[3] + LanBreWig_parameters[0]*TMath::Landau(x,LanBreWig_parameters[1],LanBreWig_parameters[2]);
    return value;
}

Double_t func_LanBreWig_fullfit(float x,Double_t *par){
    Double_t value = par[0]*TMath::BreitWigner(x,par[1],par[2]) + par[6] + par[3]*TMath::Landau(x,par[4],par[5]);
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

static float raw_fcn_BreitWigner_min;
static float raw_fcn_BreitWigner_max;
void raw_fcn_BreitWigner(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (raw_projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(raw_projection->GetBinCenter(i), raw_fcn_BreitWigner_min, raw_fcn_BreitWigner_max)*(raw_projection->GetBinContent(i)-func_BreitWigner(raw_projection->GetBinCenter(i),par))/raw_projection->GetBinError(i);
        //cout << "delta " << i << ": " << delta << endl;
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

static float fcn_Landau_K_gap_min;
static float fcn_Landau_K_gap_max;
static float fcn_Landau_rho_gap_min;
static float fcn_Landau_rho_gap_max;
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
        delta  = (enforce_interval(projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(projection->GetBinCenter(i), fcn_Landau_K_gap_max, fcn_Landau_rho_gap_min)+enforce_interval(projection->GetBinCenter(i), fcn_Landau_rho_gap_max, 1400))*(projection->GetBinContent(i)-func_Landau(projection->GetBinCenter(i),par))/(projection->GetBinError(i));
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void raw_fcn_Landau(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (raw_projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = (enforce_interval(raw_projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(raw_projection->GetBinCenter(i), fcn_Landau_K_gap_max, fcn_Landau_rho_gap_min)+enforce_interval(raw_projection->GetBinCenter(i), fcn_Landau_rho_gap_max, 1400))*(raw_projection->GetBinContent(i)-func_Landau(raw_projection->GetBinCenter(i),par))/(raw_projection->GetBinError(i));
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

static float fcn_LanBreWig_min;
static float fcn_LanBreWig_max;
void fcn_LanBreWig(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(projection->GetBinCenter(i), fcn_LanBreWig_min, fcn_LanBreWig_max)*(projection->GetBinContent(i)-func_LanBreWig(projection->GetBinCenter(i),par))/(projection->GetBinError(i));
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void raw_fcn_LanBreWig(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (raw_projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = enforce_interval(raw_projection->GetBinCenter(i), fcn_LanBreWig_min, fcn_LanBreWig_max)*(raw_projection->GetBinContent(i)-func_LanBreWig(raw_projection->GetBinCenter(i),par))/(raw_projection->GetBinError(i));
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void fcn_LanBreWig_fullfit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = (enforce_interval(projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(projection->GetBinCenter(i), fcn_Landau_K_gap_max, 1400))*(projection->GetBinContent(i)-func_LanBreWig_fullfit(projection->GetBinCenter(i),par))/(projection->GetBinError(i));
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}

void raw_fcn_LanBreWig_fullfit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (raw_projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = (enforce_interval(raw_projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(raw_projection->GetBinCenter(i), fcn_Landau_K_gap_max, 1400))*(raw_projection->GetBinContent(i)-func_LanBreWig_fullfit(raw_projection->GetBinCenter(i),par))/(raw_projection->GetBinError(i));
        chisq += delta*delta;
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl << endl;
}


string rounded(float value, int decimals) {
    float multiplier = pow(10, decimals);
    string rounded_value = to_string(round(multiplier*value)/multiplier);
    string output = "";
    bool crossed_to_decimals = false;
    if (decimals >= 0) {
        for (int i=0;i<rounded_value.length()-decimals;i++) {
            cout << output << endl;
            char symbol = rounded_value[i];
            if (symbol == '.') {
                crossed_to_decimals = true;
                if (decimals == 0) {
                    break;
                }
            } else if (crossed_to_decimals) {
                if (decimals == 0) {
                    break;
                }
                decimals -= 1;
            }
            output += symbol;
        }
    } else {
        for (int i=0;i<rounded_value.length()-decimals;i++) {
            output += rounded_value[i];
            if (rounded_value[i-decimals] == '.') {
                break;
            }
        }
    }
    return output;
}

void Monte_Carlo_4_track_glueball_analysis() {

    const string filename = "./ntuples/rho.root"; //"TOTEM43.root", kpkm.roo, 110000.root, rho.root, MinBias.root
    bool monte_carlo = false;
    if (filename == "rho.root") {
        monte_carlo = true;
    }

    const float pion_mass = 139.57039;//139.57039
    const float static_rho_mass = 774.973145; //720, 770, 743, 775.02
    float dynamic_rho_mass = static_rho_mass;
    const float static_K_mass = 498.170974;
    float dynamic_K_mass = static_K_mass;

    const float allowed_px_difference = 200;
    const float allowed_py_difference = 200;
    const float allowed_dxy_variance = 0.15;
    const float allowed_dz_variance = 0.2;
    float allowed_rho_mass_difference = 212.210888;
    const float allowed_rho_mass_difference_supercut = 50;
    
    float raw_K_radius = 37.935224;
    float K_radius = 0.65*raw_K_radius;
    const float allowed_greatest_dxy = 0.2; //0.2
    const float allowed_greatest_dz = 0.5; //0.6

    const float rho_gap_width = 212.210884;

    map<string, Double_t> results;

    TFile *file = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)file->Get("tree");
    TTreeReader Reader(tree);
    //TTreeReader Reader("tree", file);

    auto rho_masses = new TH2F("rho_masses", "Masses of potential rho particles with impact parameter cuts;m1/MeV;m2/MeV",200,200,1400,200,200,1400);
    auto rho_masses_raw = new TH2F("rho_masses_raw", ";particle 1 mass (MeV);particle 2 mass (MeV)",200,200,1400,200,200,1400);

    auto dxy_maximum_vs_rho_m1 = new TH2F("dxy_maximum_vs_rho_m1", ";dxy maximum (cm);particle 1 mass (MeV)",200,0,1,200,200,1400);
    auto dz_maximum_vs_rho_m1 = new TH2F("dz_maximum_vs_rho_m1", ";dz maximum (cm);particle 1 mass (MeV)",200,0,2,200,200,1400);

    auto dxy_maximum_rho_distribution = new TH1F("dxy_maximum_rho_distribution", "dxy maximum distribution near rho mass;dxy maximum",200,0,1);
    auto dz_maximum_rho_distribution = new TH1F("dz_maximum_rho_distribution", "dz maximum distribution near rho mass;dz maximum",200,0.2,2);

    auto raw_origin_mass = new TH1F("raw_origin_mass", "Mass of the four track origin without cuts;mass/MeV",200,500,3000);
    auto origin_mass = new TH1F("origin_mass", "Mass of the four track origin with cuts;mass/MeV",200,500,3000);
    //auto origin_mass = new TH1F("origin_mass", "Mass of the four track origin with cuts;mass/MeV",200,2000,2500);

    auto raw_origin_m_vs_rho1_m = new TH2F("raw_origin_m_vs_rho1_m", "Mass of the four track origin compared with mass of first rho particle without cuts;origin mass/MeV;rho1 mass/MeV",200,500,3000,200,200,1400);
    auto origin_m_vs_rho1_m = new TH2F("origin_m_vs_rho1_m", "Mass of the four track origin compared with mass of first rho particle with cuts;origin mass/MeV;rho1 mass/MeV",200,500,3000,200,200,1400);

    auto origin_m_vs_rho1_m_supercut = new TH2F("origin_m_vs_rho1_m_supercut", "Mass of the four track origin compared with mass of first rho particle with supercuts;origin mass/MeV;rho1 mass/MeV",200,1000,3000,200,200,1400);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kCividis);
    rho_masses->Sumw2();
    rho_masses_raw->Sumw2();
    origin_m_vs_rho1_m_supercut->Sumw2();

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

        current_event.calculate_momentum_of_diffractive_system();
        current_event.calculate_proton_momentums();

        Particle* rhos[2][2];
        current_event.reconstruct_2_from_4(rhos, 2);
        

        float raw_masses[2];
        for (int i=0; i<2; ++i) {
            for (int j=0; j<2; ++j) {
                Particle* rho = rhos[i][j];
                raw_masses[j] = rho->m;
            }
            rho_masses_raw->Fill(raw_masses[0], raw_masses[1]);

            dxy_maximum_vs_rho_m1->Fill(current_event.get_greatest_dxy(), raw_masses[1]);
            dz_maximum_vs_rho_m1->Fill(current_event.get_greatest_dz(), raw_masses[1]);
            
        }

        

        if (static_rho_mass - allowed_rho_mass_difference < rhos[0][0]->m && rhos[0][0]->m < static_rho_mass + allowed_rho_mass_difference && static_rho_mass - allowed_rho_mass_difference < rhos[0][1]->m && rhos[0][1]->m < static_rho_mass + allowed_rho_mass_difference) {
            dxy_maximum_rho_distribution->Fill(current_event.get_greatest_dxy());
            dz_maximum_rho_distribution->Fill(current_event.get_greatest_dz());
        } else if (static_rho_mass - allowed_rho_mass_difference < rhos[1][0]->m && rhos[1][0]->m < static_rho_mass + allowed_rho_mass_difference && static_rho_mass - allowed_rho_mass_difference < rhos[1][1]->m && rhos[1][1]->m < static_rho_mass + allowed_rho_mass_difference) {
            dxy_maximum_rho_distribution->Fill(current_event.get_greatest_dxy());
            dz_maximum_rho_distribution->Fill(current_event.get_greatest_dz());
        }


        if (current_event.get_greatest_dxy()>allowed_greatest_dxy) {
            continue;
        }

        if (current_event.get_greatest_dz()>allowed_greatest_dz) {
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


    auto raw_secondary_projection_canvas = new TCanvas("raw_secondary_projection_canvas","raw_secondary_projection_canvas");
    auto raw_secondary_projection = rho_masses->ProjectionY("Y_projection_raw");
    raw_secondary_projection->Draw();

    float secondary_rho_peak = 774.973145;
    float secondary_rho_min = secondary_rho_peak - 100;
    float secondary_rho_max = secondary_rho_peak + 100;

    TF1 secondary_gaussian_fit("secondary_gaussian_fit", "[0]*TMath::Gaus(x,[1],[2])", secondary_rho_min, secondary_rho_max);
    secondary_gaussian_fit.SetParameters(6.54848e+02, secondary_rho_peak, 7.43751e+01);
    raw_secondary_projection->Fit(&secondary_gaussian_fit, "","",secondary_rho_min, secondary_rho_max);
    secondary_gaussian_fit.DrawCopy("Same");

    float secondary_gaussian_fit_mean = secondary_gaussian_fit.GetParameter(1);
    float secondary_gaussian_fit_std_deviation = secondary_gaussian_fit.GetParameter(2);

    allowed_rho_mass_difference = 3*secondary_gaussian_fit_std_deviation;


    results["secondary_rho_gaussian_fit_peak"] = secondary_gaussian_fit_mean;
    results["secondary_rho_gaussian_fit_std_dev"] = secondary_gaussian_fit_std_deviation;


    TLine secondary_gaussian_line_1 = TLine(secondary_gaussian_fit_mean-3*secondary_gaussian_fit_std_deviation, 0, secondary_gaussian_fit_mean-3*secondary_gaussian_fit_std_deviation, 1.05*raw_secondary_projection->GetMaximum());
    secondary_gaussian_line_1.DrawClone();

    TLine secondary_gaussian_line_2 = TLine(secondary_gaussian_fit_mean+3*secondary_gaussian_fit_std_deviation, 0, secondary_gaussian_fit_mean+3*secondary_gaussian_fit_std_deviation, 1.05*raw_secondary_projection->GetMaximum());
    secondary_gaussian_line_2.DrawClone();


    auto dxy_maximum_distribution = new TCanvas("dxy_maximum_distribution_canvas","dxy_maximum_distribution_canvas");
    dxy_maximum_rho_distribution->Draw();

    auto dz_maximum_distribution = new TCanvas("dz_maximum_distribution_canvas","dz_maximum_distribution_canvas");
    dz_maximum_rho_distribution->Draw();


    auto raw_projections = new TCanvas("Canvas3","Canvas3");
    raw_projection = rho_masses_raw->ProjectionX("X_projection_raw");
    raw_projection->Draw();

    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    Double_t val1,err1,val2,err2,val3,err3,val4,err4,val5,err5,val6,err6,val7,err7;
    Double_t arglist[10];
    Int_t ierflg = 0;
    static Double_t vstart[7];
    static Double_t step[7];

///////////////////////////////////
//          RAW GAUSS FIT
//////////////////////////////////

    float difference = 1.;
    int i = 0;
    int repeats = 0;
    float peak1 = -1;
    float peak2 = -1;
    while (difference > i*0.0005) {

        raw_fcn_gaussian_min = dynamic_rho_mass - results["secondary_rho_gaussian_fit_std_dev"];
        raw_fcn_gaussian_max = dynamic_rho_mass + results["secondary_rho_gaussian_fit_std_dev"];

        TMinuit *gMinuit0 = new TMinuit(3);
        //gMinuit0->Command("SET PRINT -1");
        //gMinuit0->Command("SET NOWarnings");
        gMinuit0->SetFCN(raw_fcn_gaussian);

        arglist[0] = 1;

        gMinuit0->mnexcm("SET ERR", arglist ,1,ierflg);

        vstart[0] = 2000;
        vstart[1] = dynamic_rho_mass;
        vstart[2] = 50;
        step[0] = 1;
        step[1] = 1;
        step[2] = 1;
        gMinuit0->mnparm(0, "a0", vstart[0], step[0], 0,0,ierflg);
        gMinuit0->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
        gMinuit0->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);

        arglist[0] = 500;
        arglist[1] = 1.;
        gMinuit0->mnexcm("MIGRAD", arglist ,2,ierflg);

        //Double_t amin,edm,errdef;
        //Int_t nvpar,nparx,icstat;
        gMinuit0->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

        //Double_t val1,err1,val2,err2,val3,err3,val4,err4;
        gMinuit0->GetParameter(0, val1, err1);
        gMinuit0->GetParameter(1, val2, err2);
        gMinuit0->GetParameter(2, val3, err3);

        difference = abs(dynamic_rho_mass - val2);

        cout << "peak position " + to_string(i) + ": " + to_string(val2) << endl;
        cout << "previous position: " + to_string(dynamic_rho_mass) << endl;
        cout << "difference: " + to_string(difference) << endl;

        cout << endl;

        cout << "peak1: " + to_string(peak1) << endl;
        cout << "peak2: " + to_string(peak2) << endl;
        cout << "repeats: " + to_string(repeats) << endl;

        cout << endl;

        if (abs(val2 - peak2)<0.0001) {
            repeats += 1;
        } else {
            repeats = 0;
        }

        if (repeats > 10) {
            dynamic_rho_mass = (peak1 + peak2)/2;
            break;
        }

        peak2 = peak1;
        peak1 = val2;

        if (dynamic_rho_mass < val2) {
            dynamic_rho_mass += difference/3;
        } else {
            dynamic_rho_mass -= difference/3;
        }
        //dynamic_rho_mass = (dynamic_rho_mass + val2)/2;

        ++i;
        if (i > 200) {
            break;
        }
    }

    results["raw_rho_Gauss_peak"] = val2;

    TF1 raw_gaussian_fit("raw_gaussian_fit", "[0]*TMath::Gaus(x,[1],[2])", dynamic_rho_mass - results["secondary_rho_gaussian_fit_std_dev"], dynamic_rho_mass + results["secondary_rho_gaussian_fit_std_dev"]);
    raw_gaussian_fit.SetParameters(val1, val2, val3);
    raw_gaussian_fit.DrawCopy("Same");


///////////////////////////////////
//         BREIT WIGNER FIT
//////////////////////////////////
/*
    raw_fcn_BreitWigner_min = dynamic_rho_mass - rho_gap_width;
    raw_fcn_BreitWigner_max = dynamic_rho_mass + rho_gap_width;

    TMinuit *gMinuit03 = new TMinuit(3);
    gMinuit03->SetFCN(raw_fcn_BreitWigner);

    arglist[0] = 1;

    gMinuit03 ->mnexcm("SET ERR", arglist ,1,ierflg);

    vstart[0] = 1.42994e+05;
    vstart[1] = dynamic_rho_mass;
    vstart[2] = 1.17107e+02;
    step[0] = 1;
    step[1] = 0.1;
    step[2] = 1;
    gMinuit03 ->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
    gMinuit03 ->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
    gMinuit03 ->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);

    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit03 ->mnexcm("MIGRAD", arglist ,2,ierflg);

    //Double_t amin,edm,errdef;
    //Int_t nvpar,nparx,icstat;
    gMinuit03 ->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    gMinuit03 ->GetParameter(0, val1, err1);
    gMinuit03 ->GetParameter(1, val2, err2);
    gMinuit03 ->GetParameter(2, val3, err3);


    TF1 raw_BreitWigner_fit("raw_BreitWigner_fit", "[0]*TMath::BreitWigner(x,[1],[2])", raw_fcn_BreitWigner_min, raw_fcn_BreitWigner_max);
    raw_BreitWigner_fit.SetParameters(val1, val2, val3);
    raw_BreitWigner_fit.SetLineColor(6);
    //raw_BreitWigner_fit.DrawCopy("Same");

*/
    auto main = new TCanvas("Canvas4","Canvas4");
    rho_masses->Draw("Colz");
 

    auto projections = new TCanvas("Canvas5","Canvas5");
    

///////////////////////////////////
//          GAUSS FIT
//////////////////////////////////

    dynamic_rho_mass = static_rho_mass;
    difference = 1.;
    i = 0;
    repeats = 0;
    peak1 = -1;
    peak2 = -1;
    while (difference > i*0.0001) {

        projection = rho_masses->ProjectionX("X_projection", (dynamic_rho_mass-allowed_rho_mass_difference-200)/6, (dynamic_rho_mass+allowed_rho_mass_difference-200)/6);
        //projection = rho_masses->ProjectionX("X_projection");

        projection->Draw();

/////////////////////////////////////

        //peak_bin = projection->GetMaximumBin();
        //peak = projection->GetBinCenter(peak_bin);

        //fcn_BreitWigner_min = dynamic_rho_mass - 80;
        //fcn_BreitWigner_max = dynamic_rho_mass + 80;

        fcn_gaussian_min = dynamic_rho_mass - results["secondary_rho_gaussian_fit_std_dev"];
        fcn_gaussian_max = dynamic_rho_mass + results["secondary_rho_gaussian_fit_std_dev"];

        TMinuit *gMinuit1 = new TMinuit(3);
        //gMinuit1->Command("SET PRINT -1");
        //gMinuit1->Command("SET NOWarnings");
        gMinuit1->SetFCN(fcn_gaussian);

        arglist[0] = 1;

        gMinuit1->mnexcm("SET ERR", arglist ,1,ierflg);

        vstart[0] = 300;
        vstart[1] = dynamic_rho_mass;
        vstart[2] = 50;
        step[0] = 1;
        step[1] = 1;
        step[2] = 1;
        gMinuit1->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
        gMinuit1->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
        gMinuit1->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);

        arglist[0] = 500;
        arglist[1] = 1.;
        gMinuit1->mnexcm("MIGRAD", arglist ,2,ierflg);

        gMinuit1->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

        
        gMinuit1->GetParameter(0, val1, err1);
        gMinuit1->GetParameter(1, val2, err2);

        difference = abs(dynamic_rho_mass - val2);

        cout << "peak position " + to_string(i) + ": " + to_string(val2) << endl;
        cout << "previous position: " + to_string(dynamic_rho_mass) << endl;
        cout << "difference: " + to_string(difference) << endl;
        cout << endl;

        cout << "peak1: " + to_string(peak1) << endl;
        cout << "peak2: " + to_string(peak2) << endl;
        cout << "repeats: " + to_string(repeats) << endl;
        cout << endl;

        if (abs(val2 - peak2)<0.0001) {
            repeats += 1;
        } else {
            repeats = 0;
        }

        if (repeats > 10) {
            dynamic_rho_mass = (peak1 + peak2)/2;
            break;
        }

        peak2 = peak1;
        peak1 = val2;

        if (dynamic_rho_mass < val2) {
            dynamic_rho_mass += difference*(1+0.2*log(i+1))/3;
        } else {
            dynamic_rho_mass -= difference*(1+0.2*log(i+1))/3;
        }

        gMinuit1->GetParameter(2, val3, err3);

        ++i;
        if (i > 1000) {
            break;
        }
    }

    results["rho_Gauss_peak"] = val2;

    TF1 gaussian_fit("gaussian_fit", "[0]*TMath::Gaus(x,[1],[2])", dynamic_rho_mass - 60, dynamic_rho_mass + 60);
    gaussian_fit.SetParameters(val1, val2, val3);
    gaussian_fit.DrawCopy("Same");

//////////////////////
//          BreitWigner
//////////////////////
/*
    fcn_BreitWigner_min = dynamic_rho_mass-rho_gap_width;
    fcn_BreitWigner_max = dynamic_rho_mass+rho_gap_width;

    TMinuit *gMinuit4 = new TMinuit(3);
    gMinuit4->SetFCN(fcn_BreitWigner);

    arglist[0] = 1;

    gMinuit4 ->mnexcm("SET ERR", arglist ,1,ierflg);

    vstart[0] = 4.24303e+04;
    vstart[1] = dynamic_rho_mass;
    vstart[2] = 1.36910e+02;
    step[0] = 1;
    step[1] = 0.1;
    step[2] = 1;
    gMinuit4 ->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
    gMinuit4 ->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
    gMinuit4 ->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);

    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit4 ->mnexcm("MIGRAD", arglist ,2,ierflg);

    //Double_t amin,edm,errdef;
    //Int_t nvpar,nparx,icstat;
    gMinuit4 ->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    gMinuit4 ->GetParameter(0, val1, err1);
    gMinuit4 ->GetParameter(1, val2, err2);
    gMinuit4 ->GetParameter(2, val3, err3);

    TF1 BreitWigner_fit("BreitWigner_fit", "[0]*TMath::BreitWigner(x,[1],[2])", fcn_BreitWigner_min, fcn_BreitWigner_max);
    BreitWigner_fit.SetParameters(val1, val2, val3);
    BreitWigner_fit.SetLineColor(6);
    //BreitWigner_fit.DrawCopy("Same");

*/
/*
    auto proper_rho_fit_canvas = new TCanvas("proper_rho_fit_canvas","proper_rho_fit_canvas");
    

    float raw_scale = 0.33;

    auto raw_projection_scaled=new TH1D(*raw_projection);
    raw_projection_scaled->Scale(raw_scale);
    raw_projection_scaled->SetTitle("");
    raw_projection_scaled->SetYTitle("Events / 6 MeV");
    raw_projection_scaled->SetXTitle("Mass (MeV)");
    raw_projection_scaled->SetLineColor(kBlue+3);
    raw_projection_scaled->SetMarkerColor(kBlue+3);
    raw_projection_scaled->SetMarkerStyle(5);
    raw_projection_scaled->Draw("Same");


    projection->Draw("Same");
    projection->SetLineColor(kBlack);

    TF1 raw_LanBreWig_truefit_scaled_1("raw_LanBreWig_fit_scaled", "[0]*TMath::BreitWigner(x,[1],[2]) + [6] + [3]*TMath::Landau(x,[4],[5])", 200, static_K_mass-raw_K_radius);
    raw_LanBreWig_truefit_scaled_1.SetParameters(raw_scale*raw_fit_params[0], raw_fit_params[1], raw_fit_params[2], raw_scale*raw_fit_params[3], raw_fit_params[4], raw_fit_params[5], raw_scale*raw_fit_params[6]);
    raw_LanBreWig_truefit_scaled_1.SetLineColor(kBlue-4);
    raw_LanBreWig_truefit_scaled_1.SetLineStyle(10);
    raw_LanBreWig_truefit_scaled_1.DrawCopy("Same");

    TF1 raw_LanBreWig_truefit_scaled_2("raw_LanBreWig_fit_scaled", "[0]*TMath::BreitWigner(x,[1],[2]) + [6] + [3]*TMath::Landau(x,[4],[5])", static_K_mass-raw_K_radius, static_K_mass+raw_K_radius);
    raw_LanBreWig_truefit_scaled_2.SetParameters(raw_scale*raw_fit_params[0], raw_fit_params[1], raw_fit_params[2], raw_scale*raw_fit_params[3], raw_fit_params[4], raw_fit_params[5], raw_scale*raw_fit_params[6]);
    raw_LanBreWig_truefit_scaled_2.SetLineColor(kBlue-4);
    raw_LanBreWig_truefit_scaled_2.SetLineStyle(2);
    raw_LanBreWig_truefit_scaled_2.DrawCopy("Same");

    TF1 raw_LanBreWig_truefit_scaled_3("raw_LanBreWig_fit_scaled", "[0]*TMath::BreitWigner(x,[1],[2]) + [6] + [3]*TMath::Landau(x,[4],[5])", static_K_mass+raw_K_radius, 1400);
    raw_LanBreWig_truefit_scaled_3.SetParameters(raw_scale*raw_fit_params[0], raw_fit_params[1], raw_fit_params[2], raw_scale*raw_fit_params[3], raw_fit_params[4], raw_fit_params[5], raw_scale*raw_fit_params[6]);
    raw_LanBreWig_truefit_scaled_3.SetLineColor(kBlue-4);
    raw_LanBreWig_truefit_scaled_3.SetLineStyle(10);
    raw_LanBreWig_truefit_scaled_3.DrawCopy("Same");

    TF1 raw_Landau_fit_scaled("raw_Landau_fit_scaled", "[3] + [0]*TMath::Landau(x,[1],[2])", 200, 1400);
    raw_Landau_fit_scaled.SetParameters(raw_scale*raw_fit_params[3], raw_fit_params[4], raw_fit_params[5], raw_scale*raw_fit_params[6]);
    raw_Landau_fit_scaled.SetLineColor(1);
    //raw_Landau_fit_scaled.DrawCopy("Same");

    TF1 BreitWigner_fit_raw_scaled("BreitWigner_fit_raw_scaled", "[0]*TMath::BreitWigner(x,[1],[2])", 200, 1400);
    BreitWigner_fit_raw_scaled.SetParameters(raw_scale*raw_fit_params[0], raw_fit_params[1], raw_fit_params[2]);
    BreitWigner_fit_raw_scaled.SetLineColor(kCyan);
    BreitWigner_fit_raw_scaled.SetLineStyle(8);
    BreitWigner_fit_raw_scaled.DrawCopy("Same");


    TF1 LanBreWig_truefit_1("LanBreWig_truefit_1", "[0]*TMath::BreitWigner(x,[1],[2]) + [6] + [3]*TMath::Landau(x,[4],[5])", 200, static_K_mass-K_radius);
    LanBreWig_truefit_1.SetParameters(val1, val2, val3, val4, val5, val6, val7);
    LanBreWig_truefit_1.SetLineColor(kRed);
    LanBreWig_truefit_1.DrawCopy("Same");

    TF1 LanBreWig_truefit_2("LanBreWig_truefit_2", "[0]*TMath::BreitWigner(x,[1],[2]) + [6] + [3]*TMath::Landau(x,[4],[5])", static_K_mass-K_radius, static_K_mass+K_radius);
    LanBreWig_truefit_2.SetParameters(val1, val2, val3, val4, val5, val6, val7);
    LanBreWig_truefit_2.SetLineColor(kRed);
    LanBreWig_truefit_2.SetLineStyle(2);
    LanBreWig_truefit_2.DrawCopy("Same");

    TF1 LanBreWig_truefit_3("LanBreWig_truefit_3", "[0]*TMath::BreitWigner(x,[1],[2]) + [6] + [3]*TMath::Landau(x,[4],[5])", static_K_mass+K_radius, 1400);
    LanBreWig_truefit_3.SetParameters(val1, val2, val3, val4, val5, val6, val7);
    LanBreWig_truefit_3.SetLineColor(kRed-4);
    LanBreWig_truefit_3.DrawCopy("Same");

    //LanBreWig_truefit_1.DrawCopy("Same");
    //Landau_fit.DrawCopy("Same");
    BreitWigner_fit.SetLineColor(kMagenta);
    BreitWigner_fit.SetLineStyle(9);
    BreitWigner_fit.DrawCopy("Same");

    TF1 raw_K_fit_scaled("raw_K_fit", "[0]*TMath::Gaus(x,[1],[2])", static_K_mass-raw_K_radius, static_K_mass+raw_K_radius);
    raw_K_fit_scaled.SetParameters(raw_scale*raw_K_fit_values[0], raw_K_fit_values[1], raw_K_fit_values[2]);
    raw_K_fit_scaled.SetLineColor(kBlue-9);
    raw_K_fit_scaled.SetLineStyle(2);
    raw_K_fit_scaled.DrawCopy("Same");


    TLegend leg(.7,.3,.9,.9,"");
    leg.SetFillColor(0);
    leg.SetTextSize(0.025);
    leg.AddEntry(raw_projection_scaled,"raw data");
    leg.AddEntry(&raw_LanBreWig_truefit_scaled_3, "#splitline{BreitWigner fit with}{#splitline{Landau background}{on raw data}}");
    leg.AddEntry(&BreitWigner_fit_raw_scaled, "#splitline{BreitWigner fit}{on raw data}");
    leg.AddEntry(&raw_K_fit_scaled, "#splitline{Gaussian fit for}{K^{0}_{s} peak on raw data}");

    leg.AddEntry(projection,"#splitline{data with}{selection criteria}", "LE");
    leg.AddEntry(&LanBreWig_truefit_3,"#splitline{BreitWigner fit with}{Landau background}");

    leg.AddEntry(&BreitWigner_fit, "BreitWigner fit");
    leg.DrawClone("Same");
    

    CMS_lumi(proper_rho_fit_canvas, 17, 11);


    auto raw = new TCanvas("Canvas2","Canvas2");
    rho_masses_raw->Draw("Colz");

    TLine rho_masses_raw_line1 = TLine(200, dynamic_rho_mass-allowed_rho_mass_difference, 1400, dynamic_rho_mass-allowed_rho_mass_difference);
    rho_masses_raw_line1.SetLineColor(2);
    rho_masses_raw_line1.DrawClone();

    TLine rho_masses_raw_line2 = TLine(200, dynamic_rho_mass+allowed_rho_mass_difference, 1400, dynamic_rho_mass+allowed_rho_mass_difference);
    rho_masses_raw_line2.SetLineColor(2);
    rho_masses_raw_line2.DrawClone();

    TArrow arrow1(850,550,750,710,0.02,"|>");
    arrow1.SetLineWidth(3);
    arrow1.SetLineColor(0);
    //arrow1.SetFillColor(0);
    arrow1.DrawClone();

    TLatex text1(860,500,"#rho");
    text1.SetTextColor(0);
    text1.DrawClone();


    TArrow arrow2(350,350,450,450,0.02,"|>");
    arrow2.SetLineWidth(3);
    arrow2.SetLineColor(0);
    arrow2.DrawClone();

    TLatex text2(290,300,"K^{0}_{s}");
    text2.SetTextColor(0);
    text2.DrawClone();

    CMS_lumi(raw, 17, 33);


    auto dxy_maximum = new TCanvas("dxy_maximum_canvas","dxy_maximum_canvas");
    dxy_maximum_vs_rho_m1->Draw("Colz");
    CMS_lumi(dxy_maximum, 17, 33);

    TLine line21 = TLine(allowed_greatest_dxy, 200, allowed_greatest_dxy, 1400);
    line21.SetLineColor(2);
    line21.DrawClone();

    auto dz_maximum = new TCanvas("dz_maximum_canvas","dz_maximum_canvas");
    dz_maximum_vs_rho_m1->Draw("Colz");
    CMS_lumi(dz_maximum, 17, 33);

    TLine line22 = TLine(allowed_greatest_dz, 200, allowed_greatest_dz, 1400);
    line22.SetLineColor(2);
    line22.DrawClone();




    Reader.Restart();

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
        
        
        for (Particle* particle : particles) {
            //cout << "momentum: " << particle.p_t << ", charge: " << particle.charge << endl;
            cout << "p_x: " << particle->p_x << ", p_y: " << particle->p_y << ", p_z: " << particle->p_z << ", Dp: " 
            << particle->p - sqrt(pow(particle->p_x, 2)+pow(particle->p_y, 2)+pow(particle->p_z, 2)) 
            << ", Dp_t: " << particle->p_t - sqrt(pow(particle->p_x, 2)+pow(particle->p_y, 2)) << ", charge: " << particle->charge << endl;
        }
        

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

        current_event.calculate_momentum_of_diffractive_system();
        current_event.calculate_proton_momentums();

        Particle* rhos[2][2];
        current_event.reconstruct_2_from_4(rhos, 2);

        Particle* glueball1 = current_event.reconstruct_1_from_2(rhos[0][0], rhos[0][1], 3);
        Particle* glueball2 = current_event.reconstruct_1_from_2(rhos[1][0], rhos[1][1], 3);

        raw_origin_mass->Fill(glueball1->m);
        raw_origin_mass->Fill(glueball2->m);

        raw_origin_m_vs_rho1_m->Fill(glueball1->m, rhos[0][0]->m);
        raw_origin_m_vs_rho1_m->Fill(glueball2->m, rhos[1][0]->m);
        
        if (allowed_px_difference < abs(abs(current_event.PR_p[0]+current_event.PL_p[0]) - abs(current_event.diff_p[0]))) {
            continue;
        }

        if (allowed_py_difference < abs(abs(current_event.PR_p[1]+current_event.PL_p[1]) - abs(current_event.diff_p[1]))) {
            continue;
        }

        if (current_event.get_greatest_dxy()>allowed_greatest_dxy) {
            continue;
        }

        if (current_event.get_greatest_dz()>allowed_greatest_dz) {
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


        if (dynamic_rho_mass - allowed_rho_mass_difference < rhos[0][0]->m && rhos[0][0]->m < dynamic_rho_mass + allowed_rho_mass_difference && dynamic_rho_mass - allowed_rho_mass_difference < rhos[0][1]->m && rhos[0][1]->m < dynamic_rho_mass + allowed_rho_mass_difference) {
            origin_mass->Fill(glueball1->m);
            origin_m_vs_rho1_m->Fill(glueball1->m, rhos[0][0]->m);
        }

        if (dynamic_rho_mass - allowed_rho_mass_difference < rhos[1][0]->m && rhos[1][0]->m < dynamic_rho_mass + allowed_rho_mass_difference && dynamic_rho_mass - allowed_rho_mass_difference < rhos[1][1]->m && rhos[1][1]->m < dynamic_rho_mass + allowed_rho_mass_difference) {
            origin_mass->Fill(glueball2->m);
            origin_m_vs_rho1_m->Fill(glueball2->m, rhos[1][0]->m);
        }


        if (dynamic_rho_mass - allowed_rho_mass_difference_supercut < rhos[0][1]->m && rhos[0][1]->m < dynamic_rho_mass + allowed_rho_mass_difference_supercut) {
            origin_m_vs_rho1_m_supercut->Fill(glueball1->m, rhos[0][0]->m);
        }

        if (dynamic_rho_mass - allowed_rho_mass_difference_supercut < rhos[1][1]->m && rhos[1][1]->m < dynamic_rho_mass + allowed_rho_mass_difference_supercut) {
            origin_m_vs_rho1_m_supercut->Fill(glueball2->m, rhos[1][0]->m);
        }

        

        //cout << endl;
        //cout << endl;

    } 


    auto raw_origin_mass_canvas = new TCanvas("Canvas10","Canvas10");
    raw_origin_mass->Draw("Colz");

    float peak;
    float peak_bin;

    peak_bin = raw_origin_mass->GetMaximumBin();
    peak = raw_origin_mass->GetBinCenter(peak_bin);
    TF1 raw_origin_mass_fit("raw_origin_mass_fit", "[0]*TMath::Gaus(x,[1],[2])", 500, 3000);
    raw_origin_mass_fit.SetParameters(2000, peak, 100);
    //raw_origin_mass->Fit(&raw_origin_mass_fit, "","",peak-60,peak+100);

    auto origin_mass_canvas = new TCanvas("Canvas11","Canvas11");
    origin_mass->Draw("Colz");

    peak_bin = origin_mass->GetMaximumBin();
    peak = origin_mass->GetBinCenter(peak_bin);
    TF1 origin_mass_fit("origin_mass_fit", "[0]*TMath::Gaus(x,[1],[2])", 500, 3000);
    origin_mass_fit.SetParameters(2000, peak, 100);
    //origin_mass->Fit(&origin_mass_fit, "","",peak-80,peak+80);
    
    auto raw_origin_m_vs_rho1_m_canvas = new TCanvas("Canvas12","Canvas12");
    raw_origin_m_vs_rho1_m->Draw("Colz");

    auto origin_mass_canvas_canvas = new TCanvas("Canvas13","Canvas13");
    origin_m_vs_rho1_m->Draw("Colz");


    auto origin_mass_canvas_canvas_supercut = new TCanvas("Canvas14","Canvas14");
    origin_m_vs_rho1_m_supercut->Draw("Colz");

    TLine line1 = TLine(1000, dynamic_rho_mass-allowed_rho_mass_difference_supercut, 3000, dynamic_rho_mass-allowed_rho_mass_difference_supercut);
    line1.SetLineColor(2);
    line1.DrawClone();

    TLine line2 = TLine(1000, dynamic_rho_mass+allowed_rho_mass_difference_supercut, 3000, dynamic_rho_mass+allowed_rho_mass_difference_supercut);
    line2.SetLineColor(2);
    line2.DrawClone();

    auto origin_projection_supercut_canvas = new TCanvas("Canvas15","Canvas15");
    
    auto origin_projection_supercut = origin_m_vs_rho1_m_supercut->ProjectionX("origin_m_vs_rho1_m_projection_supercut", (dynamic_rho_mass-allowed_rho_mass_difference_supercut-200)/6, (dynamic_rho_mass+allowed_rho_mass_difference_supercut-200)/6);
    origin_projection_supercut->Draw();

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

    TF1 gaussian_fit("gaussian_fit", "[0]*TMath::Gaus(x,[1],[2])", 200, 1400);
    gaussian_fit.SetParameters(val1, val2, val3);

    gaussian_fit.DrawCopy("Same");
*/
}

/*
variance cuts:
724.594421 ± 2.811280
729.529297 ± 5.155312


120.558868 ± 13.010889
118.874001 ± 23.329580
*/