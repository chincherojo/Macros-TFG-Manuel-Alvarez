// fit_peak_gaussC.C
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFitResultPtr.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>
#include <memory>
#include <iostream>

// Devuelve el primer TTree del archivo, o nullptr
static TTree* getFirstTree(TFile* f) {
    if (!f) return nullptr;
    TIter it(f->GetListOfKeys());
    while (TKey* k = static_cast<TKey*>(it())) {
        if (strcmp(k->GetClassName(), "TTree") == 0)
            return static_cast<TTree*>(k->ReadObj());
    }
    return nullptr;
}

// Ajuste gaussiana + constante en rango de canales [chmin,chmax] (Channel==0).
// Devuelve true si hay éxito, y escribe mu,sigma y sus errores.
bool fit_peak_gaussC(int hour, double chmin, double chmax,
                     double &mu, double &sigma,
                     double &emu, double &esigma,
                     int nbins = 200, bool draw = true)
{
    mu = sigma = emu = esigma = std::numeric_limits<double>::quiet_NaN();
    if (hour < 0 || chmax <= chmin || nbins <= 0) {
        std::cerr << "ERROR: parametros invalidos\n"; return false;
    }

    TString fname = Form("hora%d.root", hour);
    std::unique_ptr<TFile> f(TFile::Open(fname, "READ"));
    if (!f || f->IsZombie()) { std::cerr << "ERROR: no se pudo abrir " << fname << "\n"; return false; }

    TTree* tree = getFirstTree(f.get());
    if (!tree) { std::cerr << "ERROR: no hay TTree en " << fname << "\n"; return false; }

    // Comprobar que existen ramas mínimas
    if (!tree->GetBranch("Energy")) { std::cerr << "ERROR: falta rama 'Energy'\n"; return false; }
    if (!tree->GetBranch("Channel")) { std::cerr << "ERROR: falta rama 'Channel'\n"; return false; }

    // Histograma en canales (suponiendo que 'Energy' almacena canales)
    TH1D h("h", "", nbins, chmin, chmax);
    h.GetXaxis()->SetTitle("Canal");
    h.GetYaxis()->SetTitle("Cuentas");

    TString cut = Form("Channel==0 && Energy>=%.3f && Energy<=%.3f", chmin, chmax);
    tree->Draw("Energy>>h", cut, "goff");
    if (h.GetEntries() <= 0) { std::cerr << "ERROR: histograma vacio\n"; return false; }

    // Parámetros iniciales
    int maxbin = h.GetMaximumBin();
    double mu0 = h.GetBinCenter(maxbin);
    double c0  = h.GetBinContent(h.FindBin(chmin));
    // Baseline como mínimo local aproximado en los 20% extremos
    c0 = std::min(c0, h.GetMinimum());
    double A0  = std::max(1.0, h.GetMaximum() - c0);
    double s0  = std::max( (chmax-chmin)/12.0, h.GetRMS()*0.6 );

    TF1 fgc("fgc", "[0] + [1]*exp(-0.5*((x-[2])/[3])^2)", chmin, chmax);
    fgc.SetParNames("Cte","Amp","Mu","Sigma");
    fgc.SetParameters(c0, A0, mu0, s0);
    fgc.SetParLimits(3, (chmax-chmin)/200.0, (chmax-chmin)); // sigma positiva razonable

    // Ajuste: silencioso, guarda resultado, usa HESSE; añade "E" si quieres MINOS
    TFitResultPtr r = h.Fit(&fgc, "Q S");

    mu     = fgc.GetParameter(2);
    sigma  = fgc.GetParameter(3);
    emu    = fgc.GetParError(2);
    esigma = fgc.GetParError(3);

    std::cout << "mu = " << mu << " +- " << emu
              << "   sigma = " << sigma << " +- " << esigma << "\n";

    if (draw) {
        TCanvas c("c","gauss+const",800,600);
        h.SetMarkerStyle(20);
        h.Draw("E");
        fgc.SetLineColor(kRed+1);
        fgc.Draw("SAME");
        TLegend leg(0.58,0.72,0.90,0.90);
        leg.AddEntry(&h, "Datos (Channel==0)", "lep");
        leg.AddEntry(&fgc, "Const + Gauss", "l");
        leg.AddEntry((TObject*)nullptr, Form("#mu = %.1f #pm %.1f", mu, emu), "");
        leg.AddEntry((TObject*)nullptr, Form("#sigma = %.1f #pm %.1f", sigma, esigma), "");
        leg.Draw();
        TString out = Form("fit_hour%d_%.0f_%.0f.png", hour, chmin, chmax);
        c.SaveAs(out);
    }

    return (r->Status()==0);
}

// Wrapper de ejemplo
void fit(int hour, double chmin, double chmax, int nbins=200) {
    double mu,sig,emu,esig;
    if (fit_peak_gaussC(hour, chmin, chmax, mu, sig, emu, esig, nbins, true)) {
        std::cout << "[OK] hour=" << hour
                  << "  mu=" << mu << " +- " << emu
                  << "  sigma=" << sig << " +- " << esig << "\n";
    } else {
        std::cout << "[FAIL]\n";
    }
}

// Funcion interactiva: pide numero de ajustes, luego cada rango, y ajusta con fit_peak_gaussC
void ajustar_varios() {
    int hour;
    std::cout << "hora: ";
    if (!(std::cin >> hour)) { std::cerr << "ERROR: entrada invalida\n"; return; }

    int n;
    std::cout << "numero de ajustes: ";
    if (!(std::cin >> n) || n <= 0) { std::cerr << "ERROR: numero invalido\n"; return; }

    const int nbins = 200;
    for (int i = 1; i <= n; ++i) {
        double chmin, chmax;
        std::cout << "rango " << i << " (chmin chmax): ";
        if (!(std::cin >> chmin >> chmax)) { std::cerr << "ERROR: entrada invalida\n"; return; }
        double mu, sig, emu, esig;
        bool ok = fit_peak_gaussC(hour, chmin, chmax, mu, sig, emu, esig, nbins, true);
        if (ok) {
            std::cout << "[OK] i=" << i
                      << "  mu=" << mu << " +- " << emu
                      << "  sigma=" << sig << " +- " << esig << "\n";
        } else {
            std::cout << "[FAIL] i=" << i << "  rango=[" << chmin << "," << chmax << "]\n";
        }
    }
}

