
class Particle {
    public:
        //unknown=0, pion=1, rho=2, glueball=3
        int type;
        float p;
        float p_t;
        int charge;
        float eta;
        float phi;
        float dEdx;
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
        const float prot_momentum = 6.5e+3;

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

//3.7, -3

const float high_k_curve_a = 4.1;
const float high_k_curve_b = -1.25;

const float low_k_curve_a = 2.4;
const float low_k_curve_b = -1.5;

bool is_kaon(float p, float dEdx) {
    if (dEdx < 4) {
        return false;
    }
    if (dEdx < low_k_curve_a/p + low_k_curve_b) {
        //cout << "FALSE: low " << low_k_curve_a/p + low_k_curve_b << " p: " << p << ", dEdx: " << dEdx << endl;
        return false;
    }
    if (high_k_curve_a/p + high_k_curve_b < dEdx) {
        //cout << "FALSE: high " << high_k_curve_a/p + high_k_curve_b << " p: " << p << ", dEdx: " << dEdx << endl;
        return false;
    }
    //cout << "TRUE" << endl;
    return true;
}

void dEdx_four_track_analysis() {

    const float pion_mass = 0.13957039;

    const string filenames[] = {
        "./ntuples/TOTEM40.root",
        "./ntuples/TOTEM41.root",
        "./ntuples/TOTEM42.root",
        "./ntuples/TOTEM43.root",
        "./ntuples/TOTEM20.root",
        "./ntuples/TOTEM21.root",
        "./ntuples/TOTEM22.root",
        "./ntuples/TOTEM23.root"
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

    const float min_dEdx = 5;
    const float max_dEdx = 6;
    const float min_p = 0.33;
    const float max_p = 0.5;

    auto rho_masses = new TH2F("rho_masses", ";m1/MeV;m2/MeV",200,0.2,1.4,200,0.2,1.4);
    auto dEdx_hist = new TH2F("dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    auto kaon_dEdx_hist = new TH2F("kaon_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);
    auto second_dEdx_hist = new TH2F("second_dEdx_hist", ";p (GeV);dE/dx",200,0,4,200,0,12);

    //gStyle->SetOptStat(0);
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
                particle->m = pion_mass;
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
            
            float total_charge = current_event.calculate_total_charge();
            if (total_charge != 0) {
                continue;
            }

            Particle* rhos[2][2];
            current_event.reconstruct_2_from_4(rhos, 2);

            float raw_masses[2];
            for (int i=0; i<2; ++i) {
                for (int j=0; j<2; ++j) {
                    Particle* rho = rhos[i][j];
                    raw_masses[j] = rho->m;
                }
                rho_masses->Fill(raw_masses[0], raw_masses[1]);
            }

            for (Particle* particle : particles){
                //cout << "hit " << particle->p << ", " << particle->dEdx << endl;
                dEdx_hist->Fill(particle->p, particle->dEdx);
                if (is_kaon(particle->p, particle->dEdx)) {
                    kaon_dEdx_hist->Fill(particle->p, particle->dEdx);
                }
            }

            Particle* positives[2];
            Particle* negatives[2];
            current_event.get_positives_and_negatives(positives, negatives);

            if (is_kaon(positives[0]->p, positives[0]->dEdx) && is_kaon(positives[1]->p, positives[1]->dEdx)) {
                second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
            }

            /*
            if (min_dEdx < positives[0]->dEdx && positives[0]->dEdx < max_dEdx && min_p < positives[0]->p && positives[0]->p < max_p && min_dEdx < positives[1]->dEdx && positives[1]->dEdx < max_dEdx && min_p < positives[1]->p && positives[1]->p < max_p) {

                second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
            } else if (min_dEdx < positives[1]->dEdx && positives[1]->dEdx < max_dEdx && min_p < positives[1]->p && positives[1]->p < max_p) {
                second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
            }*/

            /*
            if (min_dEdx < positives[0]->dEdx < max_dEdx && min_p < positives[0]->p < max_p) {
                second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
            } else if (min_dEdx < positives[1]->dEdx < max_dEdx && min_p < positives[1]->p < max_p) {
                second_dEdx_hist->Fill(negatives[0]->p, negatives[0]->dEdx);
                second_dEdx_hist->Fill(negatives[1]->p, negatives[1]->dEdx);
            }

            if (min_dEdx < negatives[0]->dEdx < max_dEdx && min_p < negatives[0]->p < max_p) {
                second_dEdx_hist->Fill(positives[0]->p, positives[0]->dEdx);
                second_dEdx_hist->Fill(positives[1]->p, positives[1]->dEdx);
            }
            if (min_dEdx < negatives[1]->dEdx < max_dEdx && min_p < negatives[1]->p < max_p) {
                second_dEdx_hist->Fill(positives[0]->p, positives[0]->dEdx);
                second_dEdx_hist->Fill(positives[1]->p, positives[1]->dEdx);
            }
            */
        }
    }

    auto rho_canvas = new TCanvas("rho_canvas","rho_canvas");
    rho_masses->Draw("Colz");

    auto base_canvas = new TCanvas("base_canvas","base_canvas");
    gPad->SetLogz();

    dEdx_hist->Draw("Colz");

    /*
    TLine up1 = TLine(min_p, max_dEdx, max_p, max_dEdx);
    up1.DrawClone();

    TLine down1 = TLine(min_p, min_dEdx, max_p, min_dEdx);
    down1.DrawClone();

    TLine left1 = TLine(min_p, min_dEdx, min_p, max_dEdx);
    left1.DrawClone();

    TLine right1 = TLine(max_p, min_dEdx, max_p, max_dEdx);
    right1.DrawClone();
    */

    TF1 high_k_curve1("high_k_curve1", "[0]/x+[1]");
    high_k_curve1.SetParameters(high_k_curve_a, high_k_curve_b);
    high_k_curve1.DrawCopy("Same");

    TF1 low_k_curve1("low_k_curve1", "[0]/x+[1]");
    low_k_curve1.SetParameters(low_k_curve_a, low_k_curve_b);
    low_k_curve1.DrawCopy("Same");

    TLine bar1 = TLine(low_k_curve_a/(4-low_k_curve_b), 4, high_k_curve_a/(4-high_k_curve_b), 4);
    bar1.DrawClone();

    
    auto K_dEdx_canvas = new TCanvas("K_dEdx_canvas","K_dEdx_canvas");
    gPad->SetLogz();

    kaon_dEdx_hist->Draw("Colz");

    TF1 high_k_curve2("high_k_curve2", "[0]/x+[1]");
    high_k_curve2.SetParameters(high_k_curve_a, high_k_curve_b);
    high_k_curve2.DrawCopy("Same");

    TF1 low_k_curve2("low_k_curve2", "[0]/x+[1]");
    low_k_curve2.SetParameters(low_k_curve_a, low_k_curve_b);
    low_k_curve2.DrawCopy("Same");

    TLine bar2 = TLine(low_k_curve_a/(4-low_k_curve_b), 4, high_k_curve_a/(4-high_k_curve_b), 4);
    bar2.DrawClone();


    auto second_dEdx_canvas = new TCanvas("second_dEdx_canvas","second_dEdx_canvas");
    gPad->SetLogz();

    second_dEdx_hist->Draw("Colz");

    TF1 high_k_curve3("high_k_curve3", "[0]/x+[1]");
    high_k_curve3.SetParameters(high_k_curve_a, high_k_curve_b);
    high_k_curve3.DrawCopy("Same");

    TF1 low_k_curve3("low_k_curve3", "[0]/x+[1]");
    low_k_curve3.SetParameters(low_k_curve_a, low_k_curve_b);
    low_k_curve3.DrawCopy("Same");

    TLine bar3 = TLine(low_k_curve_a/(4-low_k_curve_b), 4, high_k_curve_a/(4-high_k_curve_b), 4);
    bar3.DrawClone();

}