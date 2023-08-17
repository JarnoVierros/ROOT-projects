
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
        //unknown=0, pion=1, rho=2, glueball=3, phi=5
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
            eta = TMath::ATan(p_z/p);
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
static float fcn_Landau_phi_gap_min;
static float fcn_Landau_phi_gap_max;
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
        delta  = (enforce_interval(projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(projection->GetBinCenter(i), fcn_Landau_K_gap_max, fcn_Landau_phi_gap_min)+enforce_interval(projection->GetBinCenter(i), fcn_Landau_phi_gap_max, 1400))*(projection->GetBinContent(i)-func_Landau(projection->GetBinCenter(i),par))/(projection->GetBinError(i));
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
        delta  = (enforce_interval(raw_projection->GetBinCenter(i), 200, fcn_Landau_K_gap_min)+enforce_interval(raw_projection->GetBinCenter(i), fcn_Landau_K_gap_max, fcn_Landau_phi_gap_min)+enforce_interval(raw_projection->GetBinCenter(i), fcn_Landau_phi_gap_max, 1400))*(raw_projection->GetBinContent(i)-func_Landau(raw_projection->GetBinCenter(i),par))/(raw_projection->GetBinError(i));
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

void four_track_phiphi_analysis() {

    const string filename = "./ntuples/phi.root"; //"TOTEM43.root", kpkm.roo, 110000.root, rho.root, MinBias.root
    bool monte_carlo = true;
    if (filename == "rho.root") {
        monte_carlo = true;
    }

    const float pion_mass = 139.57039; //139.57039;//139.57039
    const float kaon_mass = 493.677;
    const float phi_mass = 1019.461;

    const float allowed_phi_mass_difference = 10;

    map<string, Double_t> results;

    TFile *file = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)file->Get("tree");
    TTreeReader Reader(tree);
    //TTreeReader Reader("tree", file);

    auto phi_masses_raw = new TH2F("phi_masses_raw", ";particle 1 mass (MeV);particle 2 mass (MeV)",200,970,1270,200,970,1270);

    auto raw_four_track_mass = new TH1F("raw_four_track_mass", ";mass (MeV)",200,1500,3500);

    auto four_track_mass = new TH1F("four_track_mass", ";mass (MeV)",200,1500,3500);


    //gStyle->SetOptStat(0);
    gStyle->SetPalette(kCividis);
    phi_masses_raw->Sumw2();
    raw_four_track_mass->Sumw2();
    four_track_mass->Sumw2();

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

        bool skip = false;

        for (int i=0;i<4;++i) {
            Particle* particle = new Particle(5);
            particle->p = 1000*trk_p[i];
            particle->p_t = 1000*trk_pt[i];
            particle->charge = trk_q[i];
            particle->eta = trk_eta[i];
            particle->phi = trk_phi[i];
            particle->m = kaon_mass;
            particle->calculate_3d_momentum();
            particle->calculate_energy();
            particles[i] = particle;

            particle->dxy = trk_dxy[i];
            particle->dz = trk_dz[i];
        }

        
        /*
        for (Particle* particle : particles) {
            cout << "momentum: " << particle->p_t << ", charge: " << particle->charge << endl;
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

        Particle* phis[2][2];
        current_event.reconstruct_2_from_4(phis, 2);

        Particle* fourtrack_1 = current_event.reconstruct_1_from_2(phis[0][0], phis[0][1], 3);
        raw_four_track_mass->Fill(fourtrack_1->m);

        Particle* fourtrack_2 = current_event.reconstruct_1_from_2(phis[1][0], phis[1][1], 3);
        raw_four_track_mass->Fill(fourtrack_2->m);

        float raw_masses[2];
        for (int i=0; i<2; ++i) {
            for (int j=0; j<2; ++j) {
                Particle* phi = phis[i][j];
                raw_masses[j] = phi->m;
            }
            phi_masses_raw->Fill(raw_masses[0], raw_masses[1]);
        }

        if (phi_mass - allowed_phi_mass_difference < phis[0][0]->m && phis[0][0]->m < phi_mass + allowed_phi_mass_difference && phi_mass - allowed_phi_mass_difference < phis[0][1]->m && phis[0][1]->m < phi_mass + allowed_phi_mass_difference) {
            four_track_mass->Fill(fourtrack_1->m);
        }

        if (phi_mass - allowed_phi_mass_difference < phis[1][0]->m && phis[1][0]->m < phi_mass + allowed_phi_mass_difference && phi_mass - allowed_phi_mass_difference < phis[1][1]->m && phis[1][1]->m < phi_mass + allowed_phi_mass_difference) {
            four_track_mass->Fill(fourtrack_2->m);
        }

        //cout << endl;
        //cout << endl;

    }


    auto raw_mass_canvas = new TCanvas("Canvas2","Canvas2");
    phi_masses_raw->Draw("Colz");

    TLine phi_line_1 = TLine(phi_mass-allowed_phi_mass_difference, phi_mass-allowed_phi_mass_difference, phi_mass+allowed_phi_mass_difference, phi_mass-allowed_phi_mass_difference);
    phi_line_1.SetLineColor(kRed);
    phi_line_1.SetLineWidth(3);
    phi_line_1.DrawClone();

    TLine phi_line_2 = TLine(phi_mass+allowed_phi_mass_difference, phi_mass-allowed_phi_mass_difference, phi_mass+allowed_phi_mass_difference, phi_mass+allowed_phi_mass_difference);
    phi_line_2.SetLineColor(kRed);
    phi_line_2.SetLineWidth(3);
    phi_line_2.DrawClone();
    
    TLine phi_line_3 = TLine(phi_mass-allowed_phi_mass_difference, phi_mass+allowed_phi_mass_difference, phi_mass+allowed_phi_mass_difference, phi_mass+allowed_phi_mass_difference);
    phi_line_3.SetLineColor(kRed);
    phi_line_3.SetLineWidth(3);
    phi_line_3.DrawClone();
    
    TLine phi_line_4 = TLine(phi_mass-allowed_phi_mass_difference, phi_mass-allowed_phi_mass_difference, phi_mass-allowed_phi_mass_difference, phi_mass+allowed_phi_mass_difference);
    phi_line_4.SetLineColor(kRed);
    phi_line_4.SetLineWidth(3);
    phi_line_4.DrawClone();


    auto raw_projections = new TCanvas("Canvas3","Canvas3");
    raw_projection = phi_masses_raw->ProjectionX("X_projection_raw");
    raw_projection->Draw();


    auto raw_four_track_mass_canvas = new TCanvas("raw_four_track_mass_canvas","raw_four_track_mass_canvas");
    raw_four_track_mass->Draw();

    auto four_track_mass_canvas = new TCanvas("four_track_mass_canvas","four_track_mass_canvas");
    four_track_mass->Draw();
}
