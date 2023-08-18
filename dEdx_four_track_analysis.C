
class Particle {
    public:
        //unknown=0, pion=1, kaons=4, rho=2, glueball=3
        int type;
        Double_t p;
        Double_t p_t;
        int charge;
        Double_t eta;
        Double_t phi;
        Double_t dEdx;
        Double_t p_x;
        Double_t p_y;
        Double_t p_z;
        Double_t m;
        Double_t E;

        Double_t dxy;
        Double_t dz;

        Particle(int type_in) {
            type = type_in;   
            m = 0;
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
            if (m == 0) {
                throw invalid_argument("mass of particle not set");
            }
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
        const Double_t prot_momentum = 6.5e+3;

        Particle** particles;
        int particle_count;

        Double_t ThxR;
        Double_t ThyR;
        Double_t ThxL;
        Double_t ThyL;

        Double_t PR_p[3] = {0,0,prot_momentum};
        Double_t PL_p[3] = {0,0,prot_momentum};


        Double_t diff_p[3] = {0,0,0}; //momentums of the diffractive system: 0=p_x, 1=p_y, 2=p_z

        Event(Particle* particles_in[], int particle_count_in) {
            particles = particles_in;
            particle_count = particle_count_in;
        }

        int calculate_total_charge() {
            Double_t total_charge = 0.;
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
            Particle* negatives[2];

            get_positives_and_negatives(positives, negatives);

            Particle* origin1_c1 = reconstruct_1_from_2(positives[0], negatives[0], type);
            origins[0][0] = origin1_c1;
            Particle* origin2_c1 = reconstruct_1_from_2(positives[1], negatives[1], type);
            origins[0][1] = origin2_c1;
            
            Particle* origin1_c2 = reconstruct_1_from_2(positives[0], negatives[1], type);
            origins[1][0] = origin1_c2;
            Particle* origin2_c2 = reconstruct_1_from_2(positives[1], negatives[0], type);
            origins[1][1] = origin2_c2;
        }

        void reconstruct_2_from_4_opposite_id(Particle* origins[2], int type) {
            if (particle_count != 4) {
                throw invalid_argument("invalid number of particles");
            }

            Particle* positives[2];
            Particle* negatives[2];

            get_positives_and_negatives(positives, negatives);

            if (positives[0]->type != negatives[0]->type && positives[1]->type != negatives[1]->type) {
                Particle* origin1 = reconstruct_1_from_2(positives[0], negatives[0], type);
                Particle* origin2 = reconstruct_1_from_2(positives[1], negatives[1], type);
                origins[0] = origin1;
                origins[1] = origin2;
            } else if (positives[0]->type != negatives[1]->type && positives[1]->type != negatives[0]->type) {
                Particle* origin1 = reconstruct_1_from_2(positives[0], negatives[1], type);
                Particle* origin2 = reconstruct_1_from_2(positives[1], negatives[0], type);
                origins[0] = origin1;
                origins[1] = origin2;
            } else {
                throw invalid_argument("The types don't match");
            }
        }

        void get_positives_and_negatives(Particle* positives[2], Particle* negatives[2]) {
            int positive_index = 0;
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

        Double_t get_total_dxy() {
            Double_t total_dxy = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_dxy += current_particle->dxy;
            }
            return total_dxy;
        }

        Double_t get_total_dz() {
            Double_t total_dz = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_dz += current_particle->dz;
            }
            return total_dz;
        }

        Double_t get_dxy_variance() {
            Double_t average_dxy = get_total_dxy()/particle_count;
            Double_t variance = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                variance += pow(average_dxy-current_particle->dxy, 2);
            }
            return variance;
        }

        Double_t get_dz_variance() {
            Double_t average_dz = get_total_dz()/particle_count;
            Double_t variance = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                variance += pow(average_dz-current_particle->dz, 2);
            }
            return variance;
        }

        Double_t get_squared_total_dxy() {
            Double_t total_squared_dxy = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_squared_dxy += pow(current_particle->dxy, 2);
            }
            return total_squared_dxy;
        }

        Double_t get_squared_total_dz() {
            Double_t total_squared_dz = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                total_squared_dz += pow(current_particle->dz, 2);
            }
            return total_squared_dz;
        }

        Double_t get_greatest_dxy() {
            Double_t maximum = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                Double_t dxy = abs(current_particle->dxy);
                if (maximum < dxy) {
                    maximum = dxy;
                }
            }
            return maximum;
        }

        Double_t get_greatest_dz() {
            Double_t maximum = 0;
            for (int i=0; i<particle_count; ++i) {
                Particle* current_particle = *(particles+i);
                Double_t dz = abs(current_particle->dz);
                if (maximum < dz) {
                    maximum = dz;
                }
            }
            return maximum;
        }
};

const Double_t high_k_curve_a = 0.336;
const Double_t high_k_curve_b = 0.14;
const Double_t high_k_curve_c = 3.04;

const Double_t low_k_curve_a = 0.064;
const Double_t low_k_curve_b = 1.33;
const Double_t low_k_curve_c = 3;

const Double_t kaon_vertical_bar = 1.2;
const Double_t kaon_horizontal_bar = 4;//3

const Double_t relaxed_kaon_vertical_bar = 1.8;

const Double_t pion_vertical_bar = 0.7;
const Double_t pion_horizontal_bar = 1.5;

const Double_t strict_pion_vertical_bar = 0.5;
const Double_t strict_pion_horizontal_bar = 3.5;

const Double_t unknown_vertical_bar = 1.5;

const Double_t two_track_mass_low_limit = 0.892-0.050;
const Double_t two_track_mass_high_limit = 0.892+0.050;

bool is_kaon(Double_t p, Double_t dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < kaon_horizontal_bar) {
        return false;
    }
    if (kaon_vertical_bar < p) {
        return false;
    }
    if (dEdx < -1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c) {
        return false;
    }
    if (-1*high_k_curve_a/(p*p)*log(high_k_curve_b*p*p) + high_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool is_pion(Double_t p, Double_t dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < pion_horizontal_bar) {
        return false;
    }
    if (pion_vertical_bar < p) {
        return false;
    }
    if (-1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool strict_is_pion(Double_t p, Double_t dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < strict_pion_horizontal_bar) {
        return false;
    }
    if (strict_pion_vertical_bar < p) {
        return false;
    }
    if (-1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool relaxed_is_kaon(Double_t p, Double_t dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < kaon_horizontal_bar) {
        return false;
    }
    if (relaxed_kaon_vertical_bar < p) {
        return false;
    }
    if (dEdx < -1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c) {
        return false;
    }
    return true;
}

bool is_unknown(Double_t p, Double_t dEdx) {
    if (dEdx < 0.1) {
        return false;
    }
    if (unknown_vertical_bar < p) {
        return true;
    }
    if (kaon_horizontal_bar < dEdx) {
        return false;
    }
    if (p < pion_vertical_bar) {
        return false;
    }
    if (-1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c < dEdx) {
        return false;
    }
    return true;
}

bool is_unknown_kaon(Double_t p, Double_t dEdx) {
    bool clear_kaon = true;
    if (dEdx < 0.1) {
        return false;
    }
    if (dEdx < kaon_horizontal_bar) {
        clear_kaon = false;
    }
    if (dEdx < -1*low_k_curve_a/(p*p)*log(low_k_curve_b*p*p) + low_k_curve_c) {
        clear_kaon = false;
    }
    if (-1*high_k_curve_a/(p*p)*log(high_k_curve_b*p*p) + high_k_curve_c < dEdx) {
        clear_kaon = false;
    }
    if (clear_kaon) {
        return true;
    }
    return is_unknown(p, dEdx);
}

bool is_unknown_pion(Double_t p, Double_t dEdx) {
    return (is_pion(p, dEdx) || is_unknown(p, dEdx));
}

Double_t distance(Double_t x, Double_t x_0, Double_t y_0, Double_t a, Double_t b, Double_t c) {
    Double_t dist = sqrt(9*pow(x-x_0, 2) + pow(-1*(a*log(b*pow(x, 2)))/pow(x, 2) + c - y_0, 2));
    //cout << "DIST: " << dist << endl;
    return dist;
}

void fcn_distance(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    f = distance(par[0], par[1], par[2], par[3], par[4], par[5]);
}

static Double_t shortest_distance_x = 0.4;
static Double_t shortest_distance_x_step = 0.1;
static TMinuit *gMinuit_distance = new TMinuit(6);
static Double_t amin_distance,edm_distance,errdef_distance;
static Int_t nvpar_distance,nparx_distance,icstat_distance;
static Double_t distance_val1,distance_err1;
static Double_t arglist_distance[10];
static Int_t ierflg_distance = 0;
Double_t shortest_distance(Double_t x, Double_t y, Double_t a, Double_t b, Double_t c){

    //cout << "parameters: " << x << ", " << y << ", " << a << ", " << b << ", " << c << endl;

    
    gMinuit_distance->Command("SET PRINT -1");
    gMinuit_distance->Command("SET NOWarnings");
    gMinuit_distance->SetFCN(fcn_distance);



    arglist_distance[0] = 1;
    gMinuit_distance->mnexcm("SET ERR",arglist_distance,1,ierflg_distance);
    gMinuit_distance->mnparm(0, "a0", shortest_distance_x, shortest_distance_x_step, 0.01,100,ierflg_distance);
    gMinuit_distance->mnparm(1, "a1", x, 0, 0,0,ierflg_distance);
    gMinuit_distance->mnparm(2, "a2", y, 0, 0,0,ierflg_distance);
    gMinuit_distance->mnparm(3, "a3", a, 0, 0,0,ierflg_distance);
    gMinuit_distance->mnparm(4, "a4", b, 0, 0,0,ierflg_distance);
    gMinuit_distance->mnparm(5, "a5", c, 0, 0,0,ierflg_distance);

    
    arglist_distance[0] = 500;
    arglist_distance[1] = 1.;
    gMinuit_distance->mnexcm("MIGRAD", arglist_distance ,2,ierflg_distance);
    gMinuit_distance->mnstat(amin_distance,edm_distance,errdef_distance,nvpar_distance,nparx_distance,icstat_distance);
    gMinuit_distance->GetParameter(0, distance_val1, distance_err1);

    shortest_distance_x = distance_val1;
    //cout << "closest x: " << distance_val1 << endl;
    Double_t dist = distance(distance_val1, x, y, a, b, c);
    //cout << "distance: " << dist << endl;
    return dist;
}

static TH2F* dEdx_data;
const Double_t max_distance = 1.0;
void fcn_curve(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

    const Int_t x_bins = 200;
    const Int_t y_bins = 200;
    Double_t chisq = 0;

    for (int i=0;i<x_bins;++i) {
        for (int j=0;j<y_bins;++j) {
            Double_t x = dEdx_data->GetXaxis()->GetBinCenter(i+1);
            Double_t y = dEdx_data->GetYaxis()->GetBinCenter(j+1);
            Int_t track_count = dEdx_data->GetBinContent(i, j);

            if (track_count == 0) {
                continue;
            }
            if (!(is_kaon(x, y))) {
                continue;
            }

            Double_t error = dEdx_data->GetBinError(i, j);
            //cout << "inputs: " << par[0] << ", " << par[1] << ", " << par[2] << endl;
            Double_t distance_from_curve = shortest_distance(x, y, par[0], par[1], par[2]);

            if (max_distance < distance_from_curve) {
                //continue;
            }
            
            Double_t delta = (pow(track_count,4)*distance_from_curve);
            chisq += delta*delta;
        }
    }
    f = chisq;
    cout << "parameters: " << par[0] << ", "<< par[1] << ", "<< par[2] << ", "<< chisq << endl;
}




void dEdx_four_track_analysis() {

    const Double_t pion_mass = 0.13957039;
    const Double_t kaon_mass = 0.493677;

    const string filenames[] = {
        "./ntuples/TOTEM40.root"/*,
        "./ntuples/TOTEM41.root",
        "./ntuples/TOTEM42.root",
        "./ntuples/TOTEM43.root",
        "./ntuples/TOTEM20.root",
        "./ntuples/TOTEM21.root",
        "./ntuples/TOTEM22.root",
        "./ntuples/TOTEM23.root"*/
    };
    /*
        "./ntuples/TOTEM40.root",
        "./ntuples/TOTEM41.root",
        "./ntuples/TOTEM42.root",
        "./ntuples/TOTEM43.root",
        "./ntuples/TOTEM20.root",
        "./ntuples/TOTEM21.root",
        "./ntuples/TOTEM22.root",
        "./ntuples/TOTEM23.root",
    */

    const Double_t min_dEdx = 5;
    const Double_t max_dEdx = 6;
    const Double_t min_p = 0.33;
    const Double_t max_p = 0.5;

    auto rho_masses = new TH2F("rho_masses", ";m1 (GeV);m2 (GeV)",200,0.2,1.4,200,0.2,1.4);
    auto two_track_mass_1 = new TH2F("two_track_mass_1", ";m1 (GeV);m2 (GeV)",100,0.6,1.6,100,0.6,1.6);
    auto four_track_mass_1 = new TH1F("four_track_mass_1", ";m (GeV)",60,1.5,3.5);
    auto four_track_mass_2 = new TH1F("four_track_mass_2", ";m (GeV)",60,1.5,3.5);

    auto dEdx_hist = new TH2F("dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    dEdx_data = dEdx_hist;
    auto kaon_dEdx_hist = new TH2F("kaon_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    auto pion_dEdx_hist = new TH2F("pion_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    auto second_dEdx_hist = new TH2F("second_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);

    four_track_mass_1->Sumw2();
    four_track_mass_2->Sumw2();

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kCividis);

    for (string filename : filenames) {
        
        cout << "Reading: " << filename << endl;

        TFile *file = TFile::Open(filename.c_str());
        TTree* tree = (TTree*)file->Get("tree");
        TTreeReader Reader(tree);

        TTreeReaderArray<Float_t> trk_dedx(Reader, "trk_dedx");
        TTreeReaderValue<Int_t> ntrk(Reader, "ntrk");

        TTreeReaderArray<Float_t> trk_p(Reader, "trk_p");
        TTreeReaderArray<Float_t> trk_pt(Reader, "trk_pt");
        TTreeReaderArray<Int_t> trk_q(Reader, "trk_q");
        TTreeReaderArray<Float_t> trk_eta(Reader, "trk_eta");
        TTreeReaderArray<Float_t> trk_phi(Reader, "trk_phi");
        while (Reader.Next()) {

            if (trk_p.GetSize() != *ntrk) {
                cout << "ALERT!" << endl;
                throw invalid_argument("ntrk is wrong!");
            }
            if (*ntrk != 4) {
                continue;
            }

            Particle* particles[4];

            for (int i=0;i<4;++i) {
                Particle* particle = new Particle(1);
                particle->p = trk_p[i];
                particle->p_t = trk_pt[i];
                particle->charge = trk_q[i];
                particle->eta = trk_eta[i];
                particle->phi = trk_phi[i];
                particle->dEdx = trk_dedx[i];
                //particle->m = pion_mass;
                particle->calculate_3d_momentum();
                //particle->calculate_energy();
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
            
            Double_t total_charge = current_event.calculate_total_charge();
            if (total_charge != 0) {
                continue;
            }

            /*
            Particle* rhos[2][2];
            current_event.reconstruct_2_from_4(rhos, 2);

            Double_t raw_masses[2];
            for (int i=0; i<2; ++i) {
                for (int j=0; j<2; ++j) {
                    Particle* rho = rhos[i][j];
                    raw_masses[j] = rho->m;
                }
                rho_masses->Fill(raw_masses[0], raw_masses[1]);
            }
            */
            for (Particle* particle : particles){
                //cout << "hit " << particle->p << ", " << particle->dEdx << endl;
                if (is_kaon(particle->p, particle->dEdx)) {
                    dEdx_hist->Fill(particle->p, particle->dEdx);
                }
                
                if (is_kaon(particle->p, particle->dEdx)) {
                    kaon_dEdx_hist->Fill(particle->p, particle->dEdx);
                }
                if (is_pion(particle->p, particle->dEdx)) {
                    pion_dEdx_hist->Fill(particle->p, particle->dEdx);
                }
            }
            
            Particle* positives[2];
            Particle* negatives[2];
            current_event.get_positives_and_negatives(positives, negatives);

            bool valid = false;
            
            if ((is_kaon(positives[0]->p, positives[0]->dEdx) && is_pion(positives[1]->p, positives[1]->dEdx)) || (is_pion(positives[0]->p, positives[0]->dEdx) && is_kaon(positives[1]->p, positives[1]->dEdx))) {
                if ((is_kaon(negatives[0]->p, negatives[0]->dEdx) && !is_kaon(negatives[1]->p, negatives[1]->dEdx)) || (is_kaon(negatives[1]->p, negatives[1]->dEdx) && !is_kaon(negatives[0]->p, negatives[0]->dEdx))) {
                    valid = true;
                }
            }

            if ((is_kaon(negatives[0]->p, negatives[0]->dEdx) && is_pion(negatives[1]->p, negatives[1]->dEdx)) || (is_pion(negatives[0]->p, negatives[0]->dEdx) && is_kaon(negatives[1]->p, negatives[1]->dEdx))) {
                if ((is_kaon(positives[0]->p, positives[0]->dEdx) && !is_kaon(positives[1]->p, positives[1]->dEdx)) || (is_kaon(positives[1]->p, positives[1]->dEdx) && !is_kaon(positives[0]->p, positives[0]->dEdx))) {
                    valid = true;
                }
            }

            if (valid) {
                int pions = 0;
                int kaons = 0;
                for (int i=0;i<4;++i) {
                    if (is_kaon(particles[i]->p, particles[i]->dEdx)) {
                        particles[i]->m = kaon_mass;
                        particles[i]->type = 4;
                        ++kaons;
                    } else {
                        particles[i]->m = pion_mass;
                        particles[i]->type = 1;
                        ++pions;
                    }
                    particles[i]->calculate_energy();
                }
                if (pions != 2 || kaons != 2) {
                    cout << "pions: " << pions << endl;
                    cout << "kaons: " << kaons << endl;
                    throw invalid_argument("Particle id went wrong!");
                }

                Particle* big_kaons[2];
                current_event.reconstruct_2_from_4_opposite_id(big_kaons, 0);

                Double_t masses[2];
                masses[0] = big_kaons[0]->m;
                masses[1] = big_kaons[1]->m;

                two_track_mass_1->Fill(masses[0], masses[1]);

                if (two_track_mass_low_limit < masses[0] && masses[0] < two_track_mass_high_limit && two_track_mass_low_limit < masses[1] && masses[1] < two_track_mass_high_limit) {
                    Particle* four_track_origin = current_event.reconstruct_1_from_2(big_kaons[0], big_kaons[1], 0);
                    four_track_mass_2->Fill(four_track_origin->m);
                }

                Particle* four_track_origin = current_event.reconstruct_1_from_2(big_kaons[0], big_kaons[1], 0);
                four_track_mass_1->Fill(four_track_origin->m);
            }

            /*
            if ((is_kaon(positives[0]->p, positives[0]->dEdx) && is_pion(positives[1]->p, positives[1]->dEdx)) || (is_pion(positives[0]->p, positives[0]->dEdx) && is_kaon(positives[1]->p, positives[1]->dEdx))) {
                second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
            }

            if ((is_kaon(negatives[0]->p, negatives[0]->dEdx) && is_pion(negatives[1]->p, negatives[1]->dEdx)) || (is_pion(negatives[0]->p, negatives[0]->dEdx) && is_kaon(negatives[1]->p, negatives[1]->dEdx))) {
                second_dEdx_hist->Fill(positives[0]->p, positives[0]->dEdx);
                second_dEdx_hist->Fill(positives[1]->p, positives[1]->dEdx);
            }
            */
            /*
            if ((is_kaon(positives[0]->p, positives[0]->dEdx) && is_pion(positives[1]->p, positives[1]->dEdx)) || (is_pion(positives[0]->p, positives[0]->dEdx) && is_kaon(positives[1]->p, positives[1]->dEdx))) {
                if (is_kaon(negatives[0]->p, negatives[0]->dEdx)) {
                    second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
                }
                if (is_kaon(negatives[1]->p, negatives[1]->dEdx)) {
                    second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                }
                
            }

            if ((is_kaon(negatives[0]->p, negatives[0]->dEdx) && is_pion(negatives[1]->p, negatives[1]->dEdx)) || (is_pion(negatives[0]->p, negatives[0]->dEdx) && is_kaon(negatives[1]->p, negatives[1]->dEdx))) {
                if (is_kaon(positives[0]->p, positives[0]->dEdx)) {
                    second_dEdx_hist->Fill(positives[1]->p, positives[1]->dEdx);
                }
                if (is_kaon(positives[1]->p, positives[1]->dEdx)) {
                    second_dEdx_hist->Fill(positives[0]->p, positives[0]->dEdx);
                }
            }
            */
           if ((is_kaon(positives[0]->p, positives[0]->dEdx) && is_pion(positives[1]->p, positives[1]->dEdx)) || (is_pion(positives[0]->p, positives[0]->dEdx) && is_kaon(positives[1]->p, positives[1]->dEdx))) {
                if (is_kaon(negatives[0]->p, negatives[0]->dEdx)) {
                    second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
                }
                if (is_kaon(negatives[1]->p, negatives[1]->dEdx)) {
                    second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                }
                
            }

            if ((is_kaon(negatives[0]->p, negatives[0]->dEdx) && is_pion(negatives[1]->p, negatives[1]->dEdx)) || (is_pion(negatives[0]->p, negatives[0]->dEdx) && is_kaon(negatives[1]->p, negatives[1]->dEdx))) {
                if (is_kaon(positives[0]->p, positives[0]->dEdx)) {
                    second_dEdx_hist->Fill(positives[1]->p, positives[1]->dEdx);
                }
                if (is_kaon(positives[1]->p, positives[1]->dEdx)) {
                    second_dEdx_hist->Fill(positives[0]->p, positives[0]->dEdx);
                }
            }
        }
    }


    
    auto bad_base_canvas = new TCanvas("bad_base_canvas","bad_base_canvas");

    dEdx_hist->Draw("Colz");

    TF1 high_k_curve1("high_k_curve1", "-[0]/(x*x)*log([1]*x*x)+[2]", 0, 1.5);
    high_k_curve1.SetParameters(high_k_curve_a, high_k_curve_b, high_k_curve_c);
    high_k_curve1.DrawCopy("Same");

    TF1 low_k_curve1("low_k_curve1", "-[0]/(x*x)*log([1]*x*x)+[2]", 0, 1.5);
    low_k_curve1.SetParameters(low_k_curve_a, low_k_curve_b, low_k_curve_c);
    low_k_curve1.DrawCopy("Same");


    auto base_canvas = new TCanvas("base_canvas","base_canvas");
    gPad->SetLogz();

    dEdx_hist->Draw("Colz");
    high_k_curve1.DrawCopy("Same");
    low_k_curve1.DrawCopy("Same");


    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    Double_t arglist[10];
    Int_t ierflg = 0;
    Double_t vstart[3];
    Double_t step[3];

    TMinuit *gMinuit = new TMinuit(3);
    //gMinuit->Command("SET PRINT -1");
    //gMinuit->Command("SET NOWarnings");
    gMinuit->SetFCN(fcn_curve);

    arglist[0] = 1;

    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    vstart[0] = 0.102268;
    vstart[1] = 0.161165;
    vstart[2] = 3.03972;
    step[0] = 0.0001;
    step[1] = 0.0001;
    step[2] = 0.0001;
    gMinuit->mnparm(0, "a0", vstart[0], step[0], 0,100,ierflg);
    gMinuit->mnparm(1, "a1", vstart[1], step[1], 0,100,ierflg);
    gMinuit->mnparm(2, "a2", vstart[2], step[2], -100,100,ierflg);

    arglist[0] = 1500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    Double_t val1,err1,val2,err2,val3,err3;
    gMinuit->GetParameter(0, val1, err1);
    gMinuit->GetParameter(1, val2, err2);
    gMinuit->GetParameter(1, val3, err3);

    cout << val1 << ", " << val2 << ", " << val3 << endl;

    TF1 fitted_k_curve1("fitted_k_curve1", "-[0]/(x*x)*log([1]*x*x)+[2]", 0, 4);
    fitted_k_curve1.SetParameters(val1, val2, val3);
    fitted_k_curve1.DrawCopy("Same");

    TF1 fitted_k_curve2("fitted_k_curve1", "-[0]/(x*x)*log([1]*x*x)+[2]", 0, 4);
    fitted_k_curve2.SetParameters(0.102268, 0.161165, 3.03972);
    fitted_k_curve2.SetLineColor(kBlue);
    fitted_k_curve2.DrawCopy("Same");

    TF1 fitted_k_curve3("fitted_k_curve1", "-[0]/(x*x)*log([1]*x*x)+[2]", 0, 4);
    fitted_k_curve1.SetParameters(0.102268, 0.161165, 3.03972);
    //fitted_k_curve1.DrawCopy("Same");

}

// kaon curve: a=0.102268, b=0.161165, c=3.03972
// pion curve: a=4.06607e-03, b=1.61135e-01, c=2.69091e+00
// higher pion curve: a=0.0793592, b=7.56146, c=2.86834