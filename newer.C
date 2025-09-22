// fit_exp_const_from_txt_hours_min.C (leyenda con chi2/ndf)
// Modelo: y(t) = B + A * exp(-t * log(2) / T12)
// Lee txt con columnas: t_s  counts  (t en segundos) y convierte a horas.
// Sin titulo, ejes en espanol, sin residuos, sin componentes por separado.

#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TString.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

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

void fit_exp_const_from_txt_hours_min(const char* txtfile,
                                      const char* out_fit_png=nullptr)
{
    std::vector<XYE> data;
    if (!read_txt_s_to_h(txtfile, data)) { std::cerr << "ERROR: sin datos\n"; return; }
    std::sort(data.begin(), data.end(), [](const XYE&a,const XYE&b){return a.x<b.x;});
    const int N = (int)data.size();
    if (N < 4) { std::cerr << "ERROR: pocos puntos\n"; return; }

    double tmin=data.front().x, tmax=data.back().x;
    double ymin=data.front().y, ymax=data.front().y;
    for (const auto& p: data){ ymin=std::min(ymin,p.y); ymax=std::max(ymax,p.y); }
    if (tmax<=tmin) tmax=tmin+1e-6;
    const double Y   = std::max(1.0, ymax);
    const double ln2 = std::log(2.0);
    const double span= std::max(1e-6, tmax - tmin);

    // inicios robustos
    int ktail=std::max(3,N/5);
    std::vector<double> tail; tail.reserve(ktail);
    for(int i=N-ktail;i<N;++i) tail.push_back(data[i].y);
    std::nth_element(tail.begin(), tail.begin()+tail.size()/2, tail.end());
    double Binit = std::max(0.0, tail[tail.size()/2]);

    int khead=std::max(3,N/3);
    double Sx=0,Sy=0,Sxx=0,Sxy=0;
    for(int i=0;i<khead;++i){
        double yy = std::max(1.0, data[i].y - Binit);
        double ly = std::log(yy);
        Sx+=data[i].x; Sy+=ly; Sxx+=data[i].x*data[i].x; Sxy+=data[i].x*ly;
    }
    double den   = khead*Sxx - Sx*Sx;
    double slope = (std::fabs(den)>0)? (khead*Sxy - Sx*Sy)/den : -1.0;
    double T12init = (slope < -1e-12)? ln2/(-slope) : std::max(span/5.0, 1e-6);
    double Ainit = std::max(1.0, ymax - Binit);

    // datos
    TGraphErrors gr(N);
    for(int i=0;i<N;++i){ gr.SetPoint(i,data[i].x,data[i].y); gr.SetPointError(i,0.0,data[i].ey); }
    gr.SetMarkerStyle(20); gr.SetMarkerSize(0.9);
    gr.SetTitle(""); // sin titulo

    // modelo y ajuste
    TF1 f("f","[0]+[1]*exp(-x*log(2)/[2])", tmin, tmax);
    f.SetParNames("B","A","T12[h]");
    f.SetParameters(Binit, Ainit, T12init);
    f.SetParLimits(0, 0.0, 10.0*Y);
    f.SetParLimits(1, 0.0, 10.0*Y);
    f.SetParLimits(2, span/100.0, span*100.0);

    TFitResultPtr r = gr.Fit(&f, "QS");
    int status = r ? r->Status() : -1;

    // grafica unica: puntos + ajuste
    gStyle->SetOptTitle(0);
    TCanvas c1("c1","", 900, 600);
    TH1F* frame = c1.DrawFrame(tmin, 0.0, tmax, std::max(1.0, ymax*1.10));
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("Tiempo (h)");
    frame->GetYaxis()->SetTitle("Cuentas");
    gr.Draw("P");
    f.SetLineWidth(2);
    f.Draw("SAME");

    // leyenda con chi2/ndf
    double chi2 = f.GetChisquare();
    int ndf = f.GetNDF();
    double chi2ndf = (ndf>0)? chi2/ndf : 0.0;

    TLegend leg(0.60, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&gr, "datos", "p");
    leg.AddEntry(&f,  "ajuste", "l");
    TString chi = TString::Format("chi2/ndf = %.3f", chi2ndf);
    leg.AddEntry((TObject*)0, chi.Data(), "");
    leg.Draw();

    TString out1 = out_fit_png ? TString(out_fit_png) : TString::Format("%s_fit_hours.png", txtfile);
    c1.SaveAs(out1);

    // salida numerica
    std::cout << "[FIT] B=" << f.GetParameter(0) << " +- " << f.GetParError(0)
              << "  A=" << f.GetParameter(1) << " +- " << f.GetParError(1)
              << "  T12[h]=" << f.GetParameter(2) << " +- " << f.GetParError(2)
              << "  chi2/ndf=" << chi2ndf
              << "  status=" << status << "\n";
}
