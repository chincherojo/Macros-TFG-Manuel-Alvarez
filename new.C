// decay_to_txt.C — exporta t[s] y cuentas en [xmin,xmax] sin histograma

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <iomanip>

static TTree* getFirstTree(TFile* f) {
    if (!f) return nullptr;
    TIter it(f->GetListOfKeys());
    while (TKey* k = static_cast<TKey*>(it()))
        if (strcmp(k->GetClassName(), "TTree") == 0)
            return static_cast<TTree*>(k->ReadObj());
    return nullptr;
}

void decay_to_txt(int firstHour, int lastHour, double xmin, double xmax)
{
    if (firstHour < 0 || lastHour < firstHour || xmax <= xmin) {
        std::cerr << "ERROR: parametros invalidos\n";
        return;
    }

    const char* PAT = "hora%d.root";   // patron de entrada
    const double DT = 3600.0;          // 1 h en segundos

    int xminI = TMath::Nint(xmin);
    int xmaxI = TMath::Nint(xmax);
    TString outname = Form("%d_%d.txt", xminI, xmaxI);

    std::ofstream out(outname.Data());
    if (!out) { std::cerr << "ERROR: no pude crear " << outname << "\n"; return; }
    out.setf(std::ios::fixed);
    out << std::setprecision(6);
    out << "# t_s\tcounts\n";

    for (int i = firstHour; i <= lastHour; ++i) {
        TString fname = Form(PAT, i);
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) { std::cerr << "WARN: no pude abrir " << fname << "\n"; continue; }

        TTree* tree = getFirstTree(f);
        if (!tree) { std::cerr << "WARN: sin TTree en " << fname << "\n"; f->Close(); continue; }

        TString cut = Form("Energy>=%.9g && Energy<=%.9g", xmin, xmax);
        if (tree->GetBranch("Channel")) cut = "Channel==0 && " + cut;

        Long64_t counts = tree->GetEntries(cut.Data());
        double t = i * DT;

        out << t << "\t" << counts << "\n";
        f->Close();
    }

    out.close();
    std::cout << "Escrito: " << outname << "\n";
}

static bool read_two_doubles(const std::string& line, double& a, double& b){
    std::string s = line;
    for(char& c: s) if(c==',' || c=='\t') c=' ';
    std::istringstream iss(s);
    return (iss >> a >> b) && (b > a);
}

void multi_decay_prompt(int firstHour, int lastHour){
    if(firstHour < 0 || lastHour < firstHour){
        std::cerr << "ERROR: rango de horas invalido\n";
        return;
    }

    std::cout << "cuantos decays vas a hacer" << std::endl;
    int N = 0;
    if(!(std::cin >> N) || N <= 0){
        std::cerr << "ERROR: N invalido\n";
        return;
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<std::pair<double,double>> ranges;
    ranges.reserve(N);

    for(int i=1;i<=N;++i){
        while(true){
            std::cout << "rango " << i << " (xmin xmax): ";
            std::cout.flush();
            std::string line;
            if(!std::getline(std::cin, line)){ std::cerr << "\nERROR: entrada interrumpida\n"; return; }
            double xmin=0, xmax=0;
            if(read_two_doubles(line, xmin, xmax)){
                ranges.emplace_back(xmin, xmax);
                break;
            }else{
                std::cerr << "entrada invalida. formato: xmin xmax, con xmax>xmin\n";
            }
        }
    }

    std::cout << "# Rangos de interes (" << ranges.size() << "):\n";
    for(size_t i=0;i<ranges.size();++i){
        std::cout << "  [" << (i+1) << "] [" << ranges[i].first << ", " << ranges[i].second << "]\n";
    }
    std::cout.flush();

    for(const auto& r: ranges){
        decay_to_txt(firstHour, lastHour, r.first, r.second);
    }
}

// fit_biexp_from_txt_hours_twofits_guard.C
// Fit1: T12_2 fijo. Fit2 (residual): T12_2 libre. Y desde 0.

#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TH1F.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#ifndef TFG_BIEXP_HELPERS
#define TFG_BIEXP_HELPERS
struct XYE { double x, y, ey; };

static inline bool read_txt_s_to_h(const char* path, std::vector<XYE>& v) {
    std::ifstream in(path);
    if (!in) return false;
    const double S2H = 1.0/3600.0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double t_s, c; if (!(ss >> t_s >> c)) continue;
        double t_h = t_s * S2H;
        double ey  = std::sqrt(std::max(0.0, c)); if (ey <= 0) ey = 1.0;
        v.push_back({t_h, c, ey});
    }
    return !v.empty();
}

static inline bool init_linear_BA1A2(double T1_h, double T2_h,
                                     const std::vector<XYE>& d, double& B, double& A1, double& A2)
{
    const double ln2 = std::log(2.0);
    double S00=0,S01=0,S02=0,S11=0,S12=0,S22=0, b0=0,b1=0,b2=0;
    for (const auto& p: d) {
        double w  = 1.0/(p.ey*p.ey);
        double e1 = std::exp(-p.x*ln2/std::max(T1_h,1e-12));
        double e2 = std::exp(-p.x*ln2/std::max(T2_h,1e-12));
        S00 += w;       S01 += w*e1;     S02 += w*e2;
        S11 += w*e1*e1; S12 += w*e1*e2;  S22 += w*e2*e2;
        b0  += w*p.y;   b1  += w*e1*p.y; b2  += w*e2*p.y;
    }
    double A[3][4] = { {S00,S01,S02,b0},{S01,S11,S12,b1},{S02,S12,S22,b2} };
    for (int c=0;c<3;++c){
        int piv=c; for(int r=c+1;r<3;++r) if (std::fabs(A[r][c])>std::fabs(A[piv][c])) piv=r;
        if (std::fabs(A[piv][c])<1e-18) return false;
        if (piv!=c) for(int k=c;k<4;++k) std::swap(A[c][k],A[piv][k]);
        double div=A[c][c]; for(int k=c;k<4;++k) A[c][k]/=div;
        for(int r=0;r<3;++r) if(r!=c){ double f=A[r][c]; for(int k=c;k<4;++k) A[r][k]-=f*A[c][k]; }
    }
    B = std::max(0.0, A[0][3]);
    A1= std::max(0.0, A[1][3]);
    A2= std::max(0.0, A[2][3]);
    return true;
}
#endif // TFG_BIEXP_HELPERS

void fit_biexp_from_txt_hours_twofits(const char* txtfile, double T12_2_fixed_h,
                                      const char* outpng1=nullptr, const char* outpng2=nullptr)
{
    if (T12_2_fixed_h <= 0) { std::cerr << "ERROR: T1/2,2 > 0\n"; return; }

    std::vector<XYE> data;
    if (!read_txt_s_to_h(txtfile, data)) { std::cerr << "ERROR: sin datos\n"; return; }
    std::sort(data.begin(), data.end(), [](const XYE&a,const XYE&b){return a.x<b.x;});
    const int N = (int)data.size(); if (N<4){ std::cerr << "ERROR: pocos puntos\n"; return; }

    double tmin=data.front().x, tmax=data.back().x;
    double ymin=data.front().y, ymax=data.front().y;
    for (const auto& p: data){ ymin=std::min(ymin,p.y); ymax=std::max(ymax,p.y); }
    if (tmax<=tmin) tmax=tmin+1e-6;
    const double Y = std::max(1.0, ymax);
    const double ln2 = std::log(2.0);
    const double span = std::max(1e-6, tmax - tmin);

    int ktail=std::max(3,N/5); std::vector<double> tail; tail.reserve(ktail);
    for(int i=N-ktail;i<N;++i) tail.push_back(data[i].y);
    std::nth_element(tail.begin(), tail.begin()+tail.size()/2, tail.end());
    double Binit = std::max(0.0, tail[tail.size()/2]);

    int khead=std::max(3,N/3);
    double Sx=0,Sy=0,Sxx=0,Sxy=0;
    for(int i=0;i<khead;++i){ double yy=std::max(1.0,data[i].y-Binit); Sx+=data[i].x; Sy+=std::log(yy); Sxx+=data[i].x*data[i].x; Sxy+=data[i].x*std::log(yy); }
    double den=khead*Sxx - Sx*Sx;
    double slope = (den!=0)? (khead*Sxy - Sx*Sy)/den : -1.0;
    double T1init = (slope<-1e-12)? ln2/(-slope) : std::max(span/5.0,1e-6);

    double A1init=0.1*Y, A2init=0.05*Y;
    init_linear_BA1A2(T1init, T12_2_fixed_h, data, Binit, A1init, A2init);

    // ===== FIT 1: total (T12_2 fijo) =====
    TGraphErrors gr(N);
    for(int i=0;i<N;++i){ gr.SetPoint(i,data[i].x,data[i].y); gr.SetPointError(i,0.0,data[i].ey); }
    gr.SetMarkerStyle(20); gr.SetMarkerSize(0.9);

    TF1 f("f","[0]+[1]*exp(-x*log(2)/[2])+[3]*exp(-x*log(2)/[4])", tmin, tmax);
    f.SetParNames("B","A1","T12_1[h]","A2","T12_2[h]");
    f.SetParameters(Binit, A1init, T1init, A2init, T12_2_fixed_h);
    f.FixParameter(4, T12_2_fixed_h);
    f.SetParLimits(0,0.0,5.0*Y);
    f.SetParLimits(1,0.0,5.0*Y);
    f.SetParLimits(3,0.0,5.0*Y);
    f.SetParLimits(2, span/100.0, span*100.0);
    gr.Fit(&f, "Q");

    TF1 fB("fB","[0]", tmin,tmax);                    fB.SetParameter(0,f.GetParameter(0));
    TF1 f1("f1","[0]*exp(-x*log(2)/[1])", tmin,tmax); f1.SetParameters(f.GetParameter(1), f.GetParameter(2));
    TF1 f2("f2","[0]*exp(-x*log(2)/[1])", tmin,tmax); f2.SetParameters(f.GetParameter(3), T12_2_fixed_h);

    f .SetLineColor(kRed+1);   f .SetLineWidth(2);
    fB.SetLineColor(kGray+2);  fB.SetLineStyle(2); fB.SetLineWidth(2);
    f1.SetLineColor(kBlue+1);  f1.SetLineStyle(3); f1.SetLineWidth(2);
    f2.SetLineColor(kGreen+2); f2.SetLineStyle(7); f2.SetLineWidth(2);

    TCanvas c1("c1","biexp fit (horas)", 900, 600);
    double ymax1 = std::max(ymax, f.Eval(tmin))*1.10;
    TH1F* frame1 = c1.DrawFrame(tmin, 0.0, tmax, std::max(1.0, ymax1));
    frame1->SetTitle(" ;Tiempo (h);Cuentas");
    gr.Draw("P");
    f.Draw("SAME"); fB.Draw("SAME"); f1.Draw("SAME"); f2.Draw("SAME");

    // Leyenda con etiquetas pedidas y chi2/ndf
    TLegend leg1(0.54,0.58,0.90,0.90);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);
    leg1.AddEntry(&gr, "Datos reales", "p");
    leg1.AddEntry(&f,  "Ajuste total", "l");
    leg1.AddEntry(&f2, "Exp. fija", "l");
    leg1.AddEntry(&f1, "Exp. libre", "l");
    leg1.AddEntry(&fB, "Constante", "l");
    leg1.AddEntry((TObject*)nullptr, Form("chi2/ndf=%.3g", f.GetChisquare()/std::max(1,f.GetNDF())), "");
    leg1.Draw();

    TString out1 = outpng1 ? TString(outpng1) : TString::Format("%s_fit_hours.png", txtfile);
    c1.SaveAs(out1);

    // ===== FIT 2: residual, T12_2 LIBRE =====
    std::vector<XYE> resid = data;
    const double A1hat = f.GetParameter(1);
    const double T1hat = f.GetParameter(2);
    const double A2hat = f.GetParameter(3);
    const double T2hat = f.GetParameter(4);
    const double B2hat = f.GetParameter(0);
    for (auto& p: resid)
        p.y = p.y - A1hat*std::exp(-p.x*ln2/std::max(T1hat,1e-12))
                    - A2hat*std::exp(-p.x*ln2/std::max(T2hat,1e-12))
                    - B2hat;

    double rmax = 0; for (auto& p: resid) rmax = std::max(rmax, p.y);
    double Y2 = std::max(1.0, rmax);

    double S00=0,S01=0,S11=0, b0=0,b1=0;
    for (const auto& p: resid) {
        double w=1.0/(p.ey*p.ey);
        double e2=std::exp(-p.x*ln2/std::max(T12_2_fixed_h,1e-12));
        S00 += w; S01 += w*e2; S11 += w*e2*e2;
        b0  += w*p.y; b1  += w*p.y*e2;
    }
    double det = S00*S11 - S01*S01;
    double B2 = std::max(0.0, (det!=0)? ( ( b0*S11 - b1*S01)/det ) : 0.0 );
    double A2i= std::max(0.0, (det!=0)? ( (-b0*S01 + b1*S00)/det ) : 0.1*Y2);
    double T2i  = std::min(std::max(T12_2_fixed_h, span/100.0), span*100.0);

    TGraphErrors gr2(N);
    for(int i=0;i<N;++i){ gr2.SetPoint(i,resid[i].x,resid[i].y); gr2.SetPointError(i,0.0,resid[i].ey); }
    gr2.SetMarkerStyle(20); gr2.SetMarkerSize(0.9);

    TF1 g("g","[0]+[1]*exp(-x*log(2)/[2])", tmin, tmax);
    g.SetParNames("B","A2","T12_2[h]");
    g.SetParameters(B2, A2i, T2i);
    g.SetParLimits(0,0.0,5.0*Y2);
    g.SetParLimits(1,0.0,5.0*Y2);
    g.SetParLimits(2, span/100.0, span*100.0);
    gr2.Fit(&g, "Q");

    TCanvas c2("c2","residual fit (horas)", 900, 600);
    double ymax2 = std::max(rmax, g.Eval(tmin))*1.10;
    TH1F* frame2 = c2.DrawFrame(tmin, 0.0, tmax, std::max(1.0, ymax2));
    frame2->SetTitle("Residual vs time;Time (h);Counts after subtracting A1-term");
    gr2.Draw("P");
    g.SetLineColor(kMagenta+1); g.SetLineWidth(2);
    g.Draw("SAME");

    TLegend leg2(0.56,0.70,0.90,0.90);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);
    leg2.AddEntry(&gr2,"Datos residuales","p");
    leg2.AddEntry(&g,"Ajuste total","l");
    leg2.AddEntry((TObject*)nullptr, Form("chi2/ndf=%.3g", g.GetChisquare()/std::max(1,g.GetNDF())), "");
    leg2.Draw();

    TString out2 = outpng2 ? TString(outpng2) : TString::Format("%s_residual_fit_hours.png", txtfile);
    c2.SaveAs(out2);

    std::cout << "[FIT1] B=" << f.GetParameter(0)
              << " A1=" << f.GetParameter(1)
              << " T12_1[h]=" << f.GetParameter(2)
              << " A2=" << f.GetParameter(3)
              << " T12_2[h](fijo)=" << T12_2_fixed_h
              << " chi2/ndf=" << f.GetChisquare()/std::max(1,f.GetNDF()) << "\n";
    std::cout << "[FIT2] B=" << g.GetParameter(0)
              << " A2=" << g.GetParameter(1)
              << " T12_2[h]=" << g.GetParameter(2)
              << " chi2/ndf=" << g.GetChisquare()/std::max(1,g.GetNDF()) << "\n";
}
