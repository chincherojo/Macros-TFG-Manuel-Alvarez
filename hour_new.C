// energy_hist_keV_from_channels.C
// Recorta por canales [xmin_ch,xmax_ch] pero dibuja X en keV.
// E_keV = (Ch - b)/a

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TString.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <memory>

static TTree* getFirstTree(TFile* f) {
    if (!f) return nullptr;
    TIter it(f->GetListOfKeys());
    while (TKey* k = (TKey*)it()) {
        if (strcmp(k->GetClassName(), "TTree")==0) return (TTree*)k->ReadObj();
    }
    return nullptr;
}

void show_hour_keV(Int_t hour,
                   Double_t xmin_ch=0, Double_t xmax_ch=4096,
                   Int_t nbins=200,
                   Int_t channel=0,
                   Double_t a=1.352, Double_t b=14.0,
                   const char* outdir="png_keV")
{
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");

    if (gSystem->AccessPathName(outdir)) gSystem->mkdir(outdir, kTRUE);

    TString fname = Form("hora%d.root", hour);
    if (gSystem->AccessPathName(fname)) { std::cerr << "[ERROR] No existe " << fname << "\n"; return; }
    std::unique_ptr<TFile> f(TFile::Open(fname));
    if (!f || f->IsZombie()) { std::cerr << "[ERROR] No pude abrir " << fname << "\n"; return; }

    TTree* tree = getFirstTree(f.get());
    if (!tree) { std::cerr << "[ERROR] Sin TTree en " << fname << "\n"; return; }

    // Convertimos los límites de canales a keV para DEFINIR el eje X
    const double xmin_keV = (xmin_ch - b)/a;
    const double xmax_keV = (xmax_ch - b)/a;

    TString hname = Form("h_keV_hora%d", hour);
    TH1D h(hname, ";Energia (keV);Cuentas", nbins, xmin_keV, xmax_keV);

    // Relleno: convertimos a keV en la expresión, pero CORTAMOS en canales
    TString drawExpr = Form("(Energy-%f)/%f>>%s", b, a, hname.Data());
    TString cutExpr  = Form("Channel==%d && Energy>=%f && Energy<=%f",
                            channel, xmin_ch, xmax_ch);

    tree->Draw(drawExpr, cutExpr, "goff");

    TCanvas c(Form("c_keV_hora%d",hour), "spectrum (keV)", 900, 600);
    c.SetTopMargin(0.12); c.SetLeftMargin(0.12);
    c.SetRightMargin(0.05); c.SetBottomMargin(0.12);

    h.SetTitle("");
    h.SetLineWidth(2);
    h.GetXaxis()->SetTitle("Energia (keV)");
    h.GetYaxis()->SetTitle("Cuentas");
    // Fuerza explícita del rango por si ROOT quiere expandir:
    h.GetXaxis()->SetRangeUser(xmin_keV, xmax_keV);

    h.Draw("HIST");

    TString out = Form("%s/hora%d_keV.png", outdir, hour);
    c.SaveAs(out);
    std::cout << "[OK] Guardado " << out << "\n";
}

// Bucle de horas [h_ini, h_fin], inclusivo. Recorta por canales y guarda en carpeta.
void batch_hours_keV(Int_t h_ini, Int_t h_fin,
                     Double_t xmin_ch=0, Double_t xmax_ch=4096,
                     Int_t nbins=200,
                     Int_t channel=0,
                     Double_t a=1.352, Double_t b=14.0,
                     const char* outdir="png_keV")
{
    if (h_ini > h_fin) std::swap(h_ini, h_fin);
    for (int h = h_ini; h <= h_fin; ++h) {
        show_hour_keV(h, xmin_ch, xmax_ch, nbins, channel, a, b, outdir);
    }
}
