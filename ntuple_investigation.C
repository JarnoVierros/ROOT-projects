
#include "CMS_lumi.C"

int get_trees(string filename, TString tree_names[]) {
    string file_path = "./ntuples/";
    if (filename == file_path+"TOTEM40.root"||filename == file_path+"TOTEM41.root"||filename == file_path+"TOTEM42.root"||filename == file_path+"TOTEM43_old.root") {
        tree_names[1] = "tree;7";
        tree_names[0] = "tree;8";
        return 1;
    } else if (filename == file_path+"TOTEM20.root"||filename == file_path+"TOTEM21.root"||filename == file_path+"TOTEM22.root"||filename == file_path+"TOTEM23.root") {
        tree_names[0] = "tree;4";
        tree_names[1] = "tree;5";
        return 2;
    } else if (filename == file_path+"TOTEM43.root") {
        tree_names[0] = "tree;1";
        tree_names[1] = "";
        return 1;
    } else if (filename == file_path+"TOTEM.root") {
        tree_names[0] = "tree;51";
        tree_names[1] = "tree;52";
        return 2;
    } else if (filename == file_path+"phi.root"||filename == file_path+"rho.root") {
        tree_names[0] = "tree;1";
        tree_names[1] = ""; 
        return 1;
    } else {
        throw invalid_argument("Unknown file");
    }

}

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

        void calculate_eta() {
            eta = TMath::ATanH(p_z/p);
        }

        void calculate_p_t() {
            p_t = sqrt(p_x*p_x + p_y*p_y);
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

static int used_bins = 0;
void fcn_LanBreWig_fullfit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    used_bins = 0;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = (enforce_interval(projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(projection->GetBinCenter(i), fcn_Landau_K_gap_max, 1400))*(projection->GetBinContent(i)-func_LanBreWig_fullfit(projection->GetBinCenter(i),par))/(projection->GetBinError(i));
        chisq += delta*delta;
        if (delta != 0) {
            ++used_bins;
        }
    }
    f = chisq;
    //cout << "chisq: " << chisq << endl;
}

static int raw_used_bins = 0;
void raw_fcn_LanBreWig_fullfit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = 200;
    Int_t i;

    raw_used_bins = 0;

    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (i=0;i<nbins; i++) {
        if (raw_projection->GetBinContent(i)<1) {
            continue;
        } 
        delta  = (enforce_interval(raw_projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(raw_projection->GetBinCenter(i), fcn_Landau_K_gap_max, 1400))*(raw_projection->GetBinContent(i)-func_LanBreWig_fullfit(raw_projection->GetBinCenter(i),par))/(raw_projection->GetBinError(i));
        chisq += delta*delta;
        if (delta != 0) {
            ++raw_used_bins;
        }
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

void ntuple_investigation() {

    const string filenames[] = {"./ntuples/TOTEM43_old.root"};
    //const string filenames[] = {"./ntuples/TOTEM40.root", "./ntuples/TOTEM41.root", "./ntuples/TOTEM42.root", "./ntuples/TOTEM43_old.root"}; //"TOTEM43.root", kpkm.roo, 110000.root, rho.root, MinBias.root
    //const string filenames[] = {"./ntuples/TOTEM20.root", "./ntuples/TOTEM21.root", "./ntuples/TOTEM22.root", "./ntuples/TOTEM23.root"};
    //{"./ntuples/TOTEM40.root", "./ntuples/TOTEM41.root", "./ntuples/TOTEM42.root", "./ntuples/TOTEM43_old.root"}
    //{"./ntuples/TOTEM20.root", "./ntuples/TOTEM21.root", "./ntuples/TOTEM22.root", "./ntuples/TOTEM23.root"}
    bool monte_carlo = false;

    const float pion_mass = 139.57039;//139.57039
    const float static_rho_mass = 724.601273; //720, 770, 743, 775.02
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

    const float four_track_rho_mass = 739;
    const float four_track_minimum_rho_mass = four_track_rho_mass - 163*(1.);
    const float four_track_maximum_rho_mass = four_track_rho_mass + 163*(1.);

    const float maximum_fourtrack_eta = 1.0; //0.65
    const float maximum_fourtrack_p_t = 5000000;

    map<string, Double_t> results;

    auto rho_masses = new TH2F("rho_masses", "Masses of potential rho particles with impact parameter cuts;m1/MeV;m2/MeV",200,200,1400,200,200,1400);
    auto rho_masses_raw = new TH2F("rho_masses_raw", ";particle 1 mass (MeV);particle 2 mass (MeV)",200,200,1400,200,200,1400);

    auto raw_four_track_mass = new TH1F("raw_four_track_mass", ";mass (MeV)",200,1000,3000);
    auto four_track_mass = new TH1F("four_track_mass", "four track mass;mass (MeV);events / 10 MeV",200,1000,3000);

    auto global_four_track_eta = new TH1F("global_four_track_eta", "four track eta of all events;events;eta",200,0,5);
    auto rho_four_track_eta = new TH1F("rho_four_track_eta", "four track eta of rho events;events;eta",200,0,5);

    auto global_four_track_pt = new TH1F("global_four_track_pt", "four track pt of all events;events;pt (MeV)",200,0,1500);
    auto rho_four_track_pt = new TH1F("rho_four_track_pt", "four track pt of rho events;events;pt (MeV)",200,0,1500);

    //gStyle->SetOptStat(0);
    gStyle->SetPalette(kCividis);
    rho_masses->Sumw2();
    rho_masses_raw->Sumw2();
    raw_four_track_mass->Sumw2();
    four_track_mass->Sumw2();
    global_four_track_eta->Sumw2();
    rho_four_track_eta->Sumw2();
    global_four_track_pt->Sumw2();
    rho_four_track_pt->Sumw2();

    int matcheds = 0;
    int unmatcheds = 0;
    int total_events = 0;
    for (string filename : filenames) {

        cout << "Reading file: " + filename << endl;

        TFile *old_file = TFile::Open(filename.c_str());

        TString trees[2];
        int tree_count = get_trees(filename, trees);
        for (int i=0; i<tree_count; ++i) {

            TTree* old_tree = (TTree*)old_file->Get(trees[i]);
            TTreeReader old_Reader(old_tree);
            //TTreeReader Reader("tree", file);

            TTreeReaderValue<UInt_t> old_Run(old_Reader, "Run");
            TTreeReaderValue<unsigned long long> old_EventNum(old_Reader, "EventNum");
            //TTreeReaderValue<UInt_t> old_LumiSection(old_Reader, "LumiSection");
            TTreeReaderArray<Int_t> old_trk_q(old_Reader, "trk_q");
            TTreeReaderValue<Float_t> old_alltrk_pt(old_Reader, "alltrk_pt");


            TFile *file = TFile::Open("./ntuples/TOTEM43.root");
            TTree* tree = (TTree*)file->Get("tree;1");
            TTreeReader Reader(tree);

            TTreeReaderValue<UInt_t> Run(Reader, "Run");
            TTreeReaderValue<unsigned long long> EventNum(Reader, "EventNum");
            //TTreeReaderValue<UInt_t> LumiSection(Reader, "LumiSection");
            TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
            TTreeReaderValue<Float_t> alltrk_pt(Reader, "alltrk_pt");
            
            int counter = 0;
            int previous_counter = 0;
            while (old_Reader.Next()) {

                if (old_trk_q.GetSize() != 4) {
                    continue;
                }
                float total_charge = 0;
                for (int i=0;i<4;++i) {
                    total_charge += old_trk_q[i];
                }
                if (total_charge != 0) {
                    continue;
                }
                ++total_events;
                /*
                TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
                TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");

                TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
                TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");

                TTreeReaderArray<Float_t> trk_dxy(Reader, "trk_dxy");
                TTreeReaderArray<Float_t> trk_dz(Reader, "trk_dz");

                TTreeReaderValue<Float_t> ThxR(Reader, "ThxR");
                TTreeReaderValue<Float_t> ThyR(Reader, "ThyR");
                TTreeReaderValue<Float_t> ThxL(Reader, "ThxL");
                TTreeReaderValue<Float_t> ThyL(Reader, "ThyL");

                TTreeReaderArray<Float_t> trk_dedx(Reader, "trk_dedx");
                */
                
                //Reader.Restart();
                counter = previous_counter;
                bool matched = false;
                counter = previous_counter;
                /*
                if (previous_counter > 0) {
                    Reader.SetEntry(previous_counter);
                    Reader.SetLocalEntry(previous_counter);
                }*/
                //Reader.SetEntriesRange(previous_counter, Reader.GetEntries());
                while (true) { //Reader.Next()

                    Reader.SetEntry(counter);
                    if (counter > previous_counter + 1000) {
                        break;
                    }
                    if (counter > Reader.GetEntries()) {
                        break;
                    }
                    //cout << Reader.GetCurrentEntry () << ", " << counter << endl;


                    if (trk_q.GetSize() != 4) {
                        ++counter;
                        continue;
                    }

                    float new_total_charge = 0;
                    for (int i=0;i<4;++i) {
                        new_total_charge += trk_q[i];
                    }
                    if (new_total_charge != 0) {
                        ++counter;
                        continue;
                    }

                    

                    if (*EventNum != *old_EventNum) {
                        ++counter;
                        continue;
                    }
                    if (*Run != *old_Run) {
                        cout << "run skip" << endl;
                        ++counter;
                        continue;
                    }

                    matched = true;
                    ++matcheds;

                    cout << counter << " (" << counter - previous_counter << "); matched " << matcheds << "/" << unmatcheds << ", event_nums: " << *EventNum << "/" << *old_EventNum << ", alltrk_pts: " << *alltrk_pt << "/" << *old_alltrk_pt << endl;
                    previous_counter = counter;
                    ++counter;
                    break;


                    /*
                    Particle* particles[4];

                    bool skip = false;

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

                        particle->dEdx = trk_dedx[i];
                    }

                    
                    bool non_pion = false;
                    for (Particle* particle : particles) {
                        if (!(is_unknown_pion((particle->p)/1000, particle->dEdx))) {
                            non_pion = true;
                        }
                    }
                    if (non_pion) {
                        continue;
                    }
                    

                    
                    for (Particle* particle : particles) {
                        cout << "momentum: " << particle->p_t << ", charge: " << particle->charge << endl;
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
                    */

                }

                if (!matched) {
                    ++unmatcheds;
                    //cout << "unmatcheds: " << unmatcheds << endl;
                }

            }
        }
    }
cout << "total events: " << total_events << endl;
}