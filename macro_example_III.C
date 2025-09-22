#include <regex>

// Input: const std::string &filename (the ROOT file that restG4 gives as output of the simulation)
// Output: this macro will display the number of "Ar37" that appeared in the event tree
void count_Ar37 (const std::string &filename){

    TRestRun *run = new TRestRun(filename);
    
    //We can print the Metadata info to check the file
    run->PrintAllMetadata();

    //We want to access the TRestGeant4Event
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    int N = 0;
    //now we can inspect all the events saved
    for (int i = 0; i < run->GetEntries(); i++){

        run->GetEntry(i);

        std::set<std::string> uniqueparticles = g4Event->GetUniqueParticles();
        
        //We print the unique particles of this event
        for (const auto& particle : uniqueparticles) {
            std::cout << particle << " ";
        }
        std::cout << std::endl;

        //We checked if some of the unique particles is "Ar37" and count it
        if (uniqueparticles.find("Ar37") != uniqueparticles.end()) {
            std::cout << "Ar37 has been produced in this event." << std::endl;
            N++;
        } else {
            std::cout << "Ar37 hasn't been produced." << std::endl;
        }
    }
    
    std::cout << "SUMMARY" << std::endl;
    std::cout << "The total number of Ar37 produced in this run is " << N << std::endl;

}

// Input: const std::string &filename (the ROOT file that restG4 gives as output of the simulation) and const std::string &isotope (the name of the isotope we want to count)
// Output: this macro will display the number of times that this isotope appeared in the event tree
void count_isotope(const std::string &filename, const std::string &isotope){

    TRestRun *run = new TRestRun(filename);

    // Imprimir la información de los metadatos para verificar el archivo
    run->PrintAllMetadata();

    // Accedemos al TRestGeant4Event
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    int totalCount = 0;

    std::map<std::string, int> isotopeCounts;
    // Construimos una expresión regular para considerar estados excitados
    // Regular expression to consider excited states, i.e Fe56[846.778]
    std::regex isotopeRegex(isotope + R"((\[\d+(\.\d+)?\])?)");

    for (int i = 0; i < run->GetEntries(); i++) {

        run->GetEntry(i);

        std::set<std::string> uniqueparticles = g4Event->GetUniqueParticles();

        //We print the unique particles of this event
        for (const auto& particle : uniqueparticles) {
            std::cout << particle << " ";
        }
        std::cout << std::endl;

        // Verify if some of these particles is the isotope that we are looking for
        for (const auto& particle : uniqueparticles) {
            if (std::regex_match(particle, isotopeRegex)) {
                isotopeCounts[particle]++;
                totalCount++;
            }
        }
    }

    std::cout << "\n========== Summary ==========\n";
    std::cout << "The total number of " << isotope << " (including excited states) produced in this run is " << totalCount << std::endl;
    for (const auto& entry : isotopeCounts) {
        std::cout << "  " << entry.first << ": " << entry.second << std::endl;
    }
    std::cout << "=========================================\n";
}


// Input: const std::string &filename (the ROOT file that restG4 gives as output of the simulation)
// Output: this macro will search for the gammas produced and make a histogram with their energy
// It will also display the name of the parent particle that has produced that emission
void plot_gammas (const std::string filename){

    TRestRun *run = new TRestRun(filename);
    
    //We can print the Metadata info to check the file
    run->PrintAllMetadata();

    //We want to access the TRestGeant4Event
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    //Here we define an arbitrary range (Think if this is the one we want)
    Double_t Emin = 100; //keV
    Double_t Emax = 3000; //keV

    //We are going to produce a ROOT histogram
    TH1D* h = new TH1D("Gammas", "Gammas emitted",500,Emin,Emax);
    //Set axis titles
    h->SetXTitle("Energy (keV)");
    h->SetYTitle("Counts");

    TString particleName;
    Double_t E_gamma;
    //now we can inspect all the events saved
    for (int i = 0; i < run->GetEntries(); i++){

        run->GetEntry(i);

        //To search the gammas, we have to access all the tracks inside this event
        for(int t = 0; t < g4Event->GetNumberOfTracks(); t++){

            auto track = g4Event->GetTrack(t);
            //we check the track "t" and get its name and initial kinetic energy
            particleName = track.GetParticleName();
            E_gamma = track.GetInitialKineticEnergy();

            //if the particle is a gamma and its energy is between the range that we are looking for, we fill the histogram
            if (particleName == "gamma" && E_gamma > Emin && E_gamma < Emax){
                h->Fill(E_gamma);

                //we can print some information
                auto parentName = (track.GetParentTrack())->GetParticleName();
                auto creatorprocess = track.GetCreatorProcess();
                std::cout << "Gamma of " << E_gamma << " keV produced via " << creatorprocess << " by " << parentName << std::endl;
                if (creatorprocess == "neutronInelastic"){
                    auto parent_track = track.GetParentTrack();
                    auto hits = parent_track->GetHits();
                    for (int n = 0; n < parent_track->GetNumberOfHits(); n++){
                        if(hits.GetHadronicTargetIsotopeName(n).size() > 0){
                            std::cout << "The hadronic process has been due to interaction with ..." << hits.GetHadronicTargetIsotopeName(n) << std::endl;
                            break;
                        }
                    }
                }
            }
        }
    }

    //To draw the histogram we have to create a canvas
    TCanvas* c = new TCanvas("c", "Gamma Energy Spectrum");
    h->Draw("same");
    c->Update();

}

//Example of usage: plot_gammas_byProcess("/home/tfg_2024_manuel/tfg/steel/steel8.root",0,10,100)
void plot_gammas_byProcess(const std::string filename, Double_t Emin, Double_t Emax, int nBins){

    TRestRun *run = new TRestRun(filename);
    run->PrintAllMetadata();
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    // Map of histograms per creator process
    std::map<std::string, TH1D*> histograms;
    //std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow};
    //int colorIndex = 0;

    gStyle->SetPalette(kBird);  // Set a palette with good contrast

    TString particleName;
    Double_t E_gamma;

    for (int i = 0; i < run->GetEntries(); i++){
        run->GetEntry(i);

        for (int t = 0; t < g4Event->GetNumberOfTracks(); t++){
            auto track = g4Event->GetTrack(t);
            particleName = track.GetParticleName();
            E_gamma = track.GetInitialKineticEnergy();

            if (particleName == "gamma" && E_gamma > Emin && E_gamma < Emax){
                TString creatorprocess = track.GetCreatorProcess();
                std::string processName = std::string(creatorprocess.Data());

                // if the process doesn't exit yet, we create a new histogram
                if (histograms.find(processName) == histograms.end()) {
                    histograms[processName] = new TH1D(Form("h_%s", creatorprocess.Data()),
                        Form("Gammas from %s", creatorprocess.Data()), nBins, Emin, Emax);
                    //histograms[processName]->SetLineColor(colors[colorIndex % colors.size()]);
                    histograms[processName]->SetXTitle("Energy (keV)");
                    histograms[processName]->SetYTitle("Counts");
                    //colorIndex++;
                }
                histograms[processName]->Fill(E_gamma);
            }
        }
    }

    // Determine the maximum count among all histograms
    Double_t maxCount = 0;

    for (const auto& [process,hist] : histograms){
        Double_t histMax = hist->GetMaximum();
        if(histMax > maxCount){
            maxCount = histMax;
        }
    }

    // plot
    TCanvas* c = new TCanvas("c", "Gamma Energy Spectrum by Process", 800, 600);
    bool first = true;
    for (auto& [process, hist] : histograms) {
        hist->GetYaxis()->SetRangeUser(0,1.1*maxCount);
        if (first) {
            hist->Draw("PLC");
            first = false;
        } else {
            hist->Draw("PLC SAME");
        }
    }
    c->BuildLegend();
    c->Update();
}

// we have to introduce the process we are interested in (RadioactiveDecay, neutronInelastic)
// we have an option to differenciate or not the excited states of the parent isotopes (if it's the case). The default value is false
// example of usage: plot_gammas_byParent("/home/tfg_2024_manuel/tfg/steel/steel8.root","RadioactiveDecay",0,10,100,false)
void plot_gammas_byParent(const std::string& filename, const std::string& process, Double_t Emin, Double_t Emax, int nBins, bool excited_states = false) {
    TRestRun *run = new TRestRun(filename);
    run->PrintAllMetadata();
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    // Map of histograms per process and parent track
    std::map<std::string, TH1D*> histograms;
    
    gStyle->SetPalette(kRainBow);  // Set a palette with good contrast

    TString particleName;
    Double_t E_gamma;

    for (int i = 0; i < run->GetEntries(); i++) {
        run->GetEntry(i);

        for (int t = 0; t < g4Event->GetNumberOfTracks(); t++) {
            auto track = g4Event->GetTrack(t);
            particleName = track.GetParticleName();
            E_gamma = track.GetInitialKineticEnergy();

            if (particleName == "gamma" && E_gamma > Emin && E_gamma < Emax) {
                TString creatorprocess = track.GetCreatorProcess();

                if (std::string(creatorprocess.Data()) != process) continue;

                std::string parentStr = "Unknown";

                if (creatorprocess == "neutronInelastic") {
                    auto parent_track = track.GetParentTrack();
                    auto hits = parent_track->GetHits();
                    for (int n = 0; n < parent_track->GetNumberOfHits(); n++) {
                        std::string isotope = hits.GetHadronicTargetIsotopeName(n);

                        if (!isotope.empty()) {
                            if (!excited_states) {
                                size_t pos = isotope.find("[");
                                if (pos != std::string::npos) {
                                    isotope = isotope.substr(0, pos);
                                }
                            }
                            parentStr = isotope;
                            break;
                        }
                    }
                } else if (creatorprocess == "RadioactiveDecay") {

                    auto parentTrack = track.GetParentTrack();
                    TString parentName = parentTrack ? parentTrack->GetParticleName() : "Unknown";
                    std::string isotope = std::string(parentName.Data());
                    size_t pos = isotope.find("[");
                    if (pos != std::string::npos) {
                        isotope = isotope.substr(0, pos);
                    }
                    parentStr = isotope;

                } else {
                    auto parentTrack = track.GetParentTrack();
                    TString parentName = parentTrack ? parentTrack->GetParticleName() : "Unknown";
                    parentStr = std::string(parentName.Data());
                }

                // Create histogram if it doesn't exist for this parent/isotope
                if (histograms.find(parentStr) == histograms.end()) {
                    histograms[parentStr] = new TH1D(Form("h_%s_%s", creatorprocess.Data(), parentStr.c_str()),
                        Form("Gammas from %s (Parent/Isotope: %s)", creatorprocess.Data(), parentStr.c_str()), nBins, Emin, Emax);
                    
                    //int color = colorIndex % nColors;  
                    //histograms[parentStr]->SetLineColor(gStyle->GetColorPalette(color));
                    //colorIndex++;                        
                    histograms[parentStr]->SetXTitle("Energy (keV)");
                    histograms[parentStr]->SetYTitle("Counts");
                }
                histograms[parentStr]->Fill(E_gamma);
            }
        }
    }

    // Determine maximum count for all histograms (to adjust the y-range)
    Double_t maxCount = 0;
    for (const auto& [parent, hist] : histograms) {
        Double_t histMax = hist->GetMaximum();
        if (histMax > maxCount) {
            maxCount = histMax;
        }
    }

    // Plot
    TCanvas* c = new TCanvas("c_gammas", "Gamma Spectrum by Parent/Isotope", 1200, 800);
    bool first = true;
    for (auto& [parent, hist] : histograms) {
        if (first) {
            hist->GetYaxis()->SetRangeUser(0, 1.1 * maxCount);
            hist->Draw("PLC");
            first = false;
        } else {
            hist->Draw("PLC SAME");
        }
    }
    c->BuildLegend();
    c->Update();

    // Summary
    std::cout << "\n========== Summary ==========\n";
    int totalGammas = 0;
    for (const auto& [parent, hist] : histograms) {
        int parentGammas = hist->GetEntries();
        totalGammas += parentGammas;
        std::cout << "Parent/Isotope: " << parent << " - Gammas: " << parentGammas << std::endl;
    }
    std::cout << "Total gammas: " << totalGammas << std::endl;
    std::cout << "=========================================\n";
}


// we have to introduce the process we are interested in (RadioactiveDecay or neutronInelastic) and the parent/target isotope
// we are not differenciating excited states of each isotope
// Example of usage: plot_gammas_Element_Process("/home/tfg_2024_manuel/tfg/steel/steel8.root","neutronInelastic","Fe",1000,2000,500)
void plot_gammas_Element_Process(const std::string& filename, const std::string& process, const std::string& element, Double_t Emin, Double_t Emax, int nBins) {
    TRestRun *run = new TRestRun(filename);
    run->PrintAllMetadata();
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    // We verify the process
    if (process != "neutronInelastic" && process != "RadioactiveDecay") {
        std::cerr << "Error: Proceso no soportado. Debe ser 'neutronInelastic' o 'RadioactiveDecay'." << std::endl;
        return;
    }

    //Map of  histograms per isotope
    std::map<std::string, TH1D*> histograms;
    TString particleName;
    Double_t E_gamma;

    gStyle->SetPalette(kRainBow);  // Set a palette with good contrast

    for (int i = 0; i < run->GetEntries(); i++) {
        run->GetEntry(i);

        for (int t = 0; t < g4Event->GetNumberOfTracks(); t++) {
            auto track = g4Event->GetTrack(t);
            particleName = track.GetParticleName();
            E_gamma = track.GetInitialKineticEnergy();

            if (particleName == "gamma" && E_gamma > Emin && E_gamma < Emax) {
                TString creatorprocess = track.GetCreatorProcess();
                if (std::string(creatorprocess.Data()) != process) continue;

                std::string parentStr = "Unknown";

                if (creatorprocess == "neutronInelastic") {
                    auto parent_track = track.GetParentTrack();
                    auto hits = parent_track->GetHits();
                    for (int n = 0; n < parent_track->GetNumberOfHits(); n++) {
                        std::string foundIsotope = hits.GetHadronicTargetIsotopeName(n);
                        if (!foundIsotope.empty()) {
                            size_t pos = foundIsotope.find("[");
                            if (pos != std::string::npos) {
                                foundIsotope = foundIsotope.substr(0, pos);
                            }
                            if (foundIsotope.find(element) == 0) {
                                parentStr = foundIsotope;
                                break;
                            }
                        }
                    }
                } else {
                    auto parentTrack = track.GetParentTrack();
                    TString parentName = parentTrack ? parentTrack->GetParticleName() : "Unknown";
                    std::string foundIsotope = std::string(parentName.Data());
                    size_t pos = foundIsotope.find("[");
                    if (pos != std::string::npos) {
                        foundIsotope = foundIsotope.substr(0, pos);
                    }
                    parentStr = foundIsotope;
                }

                if (parentStr.find(element) == 0) {
                    if (histograms.find(parentStr) == histograms.end()) {
                        histograms[parentStr] = new TH1D(Form("h_%s_%s", process.c_str(), parentStr.c_str()),
                            Form("Gamma Spectrum from %s (Isotope: %s)", process.c_str(), parentStr.c_str()), nBins, Emin, Emax);
                        histograms[parentStr]->SetLineColor(histograms.size() + 1);
                        histograms[parentStr]->SetXTitle("Energy (keV)");
                        histograms[parentStr]->SetYTitle("Counts");
                    }
                    histograms[parentStr]->Fill(E_gamma);
                }
            }
        }
    }

    Double_t maxCount = 0;
    for (const auto& [isotope, hist] : histograms) {
        Double_t histMax = hist->GetMaximum();
        if (histMax > maxCount) {
            maxCount = histMax;
        }
    }

    TCanvas* c = new TCanvas("c_gammas", "Gamma Spectrum by Element and Isotope", 1200, 800);
    bool first = true;
    for (auto& [isotope, hist] : histograms) {
        if (first) {
            hist->GetYaxis()->SetRangeUser(0, 1.1 * maxCount);
            hist->Draw("PLC");
            first = false;
        } else {
            hist->Draw("PLC SAME");
        }
    }
    c->BuildLegend();
    c->Update();

    std::cout << "\n========== Summary ==========\n";
    int totalGammas = 0;
    for (const auto& [isotope, hist] : histograms) {
        int isotopeGammas = hist->GetEntries();
        totalGammas += isotopeGammas;
        std::cout << "Isotope: " << isotope << " - Gammas: " << isotopeGammas << std::endl;
    }
    std::cout << "Total gammas: " << totalGammas << std::endl;
    std::cout << "=========================================" << std::endl;
}


// we have to introduce the process we are interested in (RadioactiveDecay or neutronInelastic) and the parent/target isotope
// we are not differenciating excited states of each isotope
void plot_gammas_Isotope_Process(const std::string& filename, const std::string& process, const std::string& isotope, Double_t Emin, Double_t Emax, int nBins) {
    TRestRun *run = new TRestRun(filename);
    run->PrintAllMetadata();
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    TH1D* histogram = new TH1D(Form("h_%s_%s", process.c_str(), isotope.c_str()),
        Form("Gamma Spectrum from %s (Isotope: %s)", process.c_str(), isotope.c_str()), nBins, Emin, Emax);
    histogram->SetLineColor(kBlue);
    histogram->SetXTitle("Energy (keV)");
    histogram->SetYTitle("Counts");

    TString particleName;
    Double_t E_gamma;

    if (process != "neutronInelastic" && process != "RadioactiveDecay") {
        std::cerr << "Error: Proceso no soportado. Debe ser 'neutronInelastic' o 'RadioactiveDecay'." << std::endl;
        return;
    }

    for (int i = 0; i < run->GetEntries(); i++) {
        run->GetEntry(i);

        for (int t = 0; t < g4Event->GetNumberOfTracks(); t++) {
            auto track = g4Event->GetTrack(t);
            particleName = track.GetParticleName();
            E_gamma = track.GetInitialKineticEnergy();

            if (particleName == "gamma" && E_gamma > Emin && E_gamma < Emax) {
                TString creatorprocess = track.GetCreatorProcess();
                if (std::string(creatorprocess.Data()) != process) continue;

                std::string parentStr = "Unknown";

                
                if (creatorprocess == "neutronInelastic") {
                    auto parent_track = track.GetParentTrack();
                    auto hits = parent_track->GetHits();
                    for (int n = 0; n < parent_track->GetNumberOfHits(); n++) {
                        std::string foundIsotope = hits.GetHadronicTargetIsotopeName(n);
                        if(!foundIsotope.empty()){
                            size_t pos = foundIsotope.find("[");
                            if (pos != std::string::npos) {
                                foundIsotope = foundIsotope.substr(0, pos); 
                            }
                            
                            if (foundIsotope == isotope) {
                                parentStr = foundIsotope;
                                break;
                            }
                        }
                    }
                } else if (creatorprocess == "RadioactiveDecay") {
                    auto parentTrack = track.GetParentTrack();
                    TString parentName = parentTrack ? parentTrack->GetParticleName() : "Unknown";
                    std::string foundIsotope = std::string(parentName.Data());
                    size_t pos = foundIsotope.find("[");
                    if (pos != std::string::npos) {
                        foundIsotope = foundIsotope.substr(0, pos); 
                    }
                    parentStr = foundIsotope;
                }

                if (parentStr.find(isotope) == 0) {
                    histogram->Fill(E_gamma);
                }
            }
        }
    }


    TCanvas* c = new TCanvas("c_gammas_isotope", "Gamma Spectrum by Process and Isotope", 800, 600);
    histogram->Draw();
    c->Update();

    std::cout << "\n========== Summary ==========" << std::endl;
    int totalGammas = histogram->GetEntries();
    std::cout << "Process: " << process << " - Isotope: " << isotope << " - Gammas: " << totalGammas << std::endl;
    std::cout << "=========================================" << std::endl;
}


void plot_gammas_new(const std::string& filename,
                                  Double_t Emin,
                                  Double_t Emax,
                                  int nBins,
                                  bool excited_states = false) {
    // Sin caja de estadísticos ni título global
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Open the run
    TRestRun* run = new TRestRun(filename);
    TRestGeant4Event* g4Event = static_cast<TRestGeant4Event*>(run->GetInputEvent());

    // Isotopes of interest and their drawing colors
    const std::vector<std::string> isotopes = {"K41", "Ca42", "Ca43", "Fe56", "Fe58", "Cl37"};
    const std::vector<int> colors          = {kRed + 1, kOrange + 7, kGreen + 2, kAzure + 1, kViolet + 5};

    std::map<std::string, TH1D*> hist;

    // Prepare one histogram per isotope
    for (size_t i = 0; i < isotopes.size(); ++i) {
        const std::string& iso = isotopes[i];
        hist[iso] = new TH1D(("h_" + iso).c_str(), "", nBins, Emin, Emax); // sin título
        hist[iso]->SetLineColor(colors[i]);
        hist[iso]->SetLineWidth(2);
        hist[iso]->SetXTitle("Energia (keV)"); // ejes en espanol sin tildes
        hist[iso]->SetYTitle("Cuentas");
        hist[iso]->SetStats(0);
    }

    // Loop over events
    for (Long64_t ievt = 0; ievt < run->GetEntries(); ++ievt) {
        run->GetEntry(ievt);

        // Loop over tracks
        for (int t = 0; t < g4Event->GetNumberOfTracks(); ++t) {
            auto track = g4Event->GetTrack(t);

            if (track.GetParticleName() != "gamma") continue;                 // Only gammas
            Double_t E = track.GetInitialKineticEnergy();
            if (E < Emin || E > Emax) continue;                                // Energy window
            if (track.GetCreatorProcess() != "RadioactiveDecay") continue;     // Only RadioactiveDecay

            // Identify parent isotope
            std::string iso = "Unknown";
            if (auto parent = track.GetParentTrack()) {
                iso = parent->GetParticleName().Data();
                size_t pos = iso.find("[");
                if (pos != std::string::npos && !excited_states) {
                    iso = iso.substr(0, pos);  // Strip excited state suffix
                }
            }

            // Fill only if this isotope is in the selected list
            if (hist.count(iso)) {
                hist[iso]->Fill(E);
            }
        }
    }

    // Determine global maximum to set Y-axis range
    Double_t maxCount = 0;
    for (auto& kv : hist) {
        maxCount = std::max(maxCount, kv.second->GetMaximum());
    }

    // Draw
    TCanvas* c = new TCanvas("c_iso", "", 1200, 800); // sin título
    bool first = true;

    // Construimos una leyenda manual (sin espacio extra a la derecha)
    TLegend* leg = new TLegend(0.68, 0.70, 0.88, 0.90); // caja mas estrecha
    leg->SetFillStyle(0);   // transparente
    leg->SetBorderSize(0);  // sin borde
    leg->SetMargin(0.12);   // menos margen interno (reduce "aire" a la derecha)
    leg->SetTextSize(0.030);

    for (const auto& iso : isotopes) {
        auto* h = hist[iso];
        if (h->GetEntries() == 0) continue;  // Skip empty spectra

        if (first) {
            h->GetYaxis()->SetRangeUser(0, 1.1 * maxCount);
            h->Draw("HIST");
            first = false;
        } else {
            h->Draw("HIST SAME");
        }
        leg->AddEntry(h, iso.c_str(), "l");
    }

    leg->Draw();
    c->Update();

    // Console summary
    std::cout << "\n==== Summary ====\n";
    int total = 0;
    for (const auto& iso : isotopes) {
        auto* h = hist[iso];
        std::cout << iso << " : " << h->GetEntries() << "\n";
        total += h->GetEntries();
    }
    std::cout << "Total gammas: " << total << "\n";
    std::cout << "==================\n";
}
