ROOT::Math::VavilovAccurate vav;
TFile *out_rootfile;

void Landau_Vavilov_comparison(){

  TF1 *vav_kappa_0p005 = new TF1("vav_kappa_0p005", "vav.Pdf(x, 0.005, 0.5)", -20., 200.);
  TF1 *vav_kappa_0p002 = new TF1("vav_kappa_0p002", "vav.Pdf(x, 0.002, 0.5)", -20., 200.);
  TF1 *vav_kappa_0p001 = new TF1("vav_kappa_0p001", "vav.Pdf(x, 0.001, 0.5)", -20., 200.);

  TF1 *landau_c_1 = new TF1("landau_c_1", "TMath::Landau(x, 0., 1.)", -20., 200.);
  TF1 *landau_c_piover2 = new TF1("landau_c_1", "TMath::Landau(x, 0., TMath::Pi() / 2.)", -20., 200.);

  double xmin = -20.;
  double xmax = 50.;
  double ymax = 0.25;

  TCanvas *c = new TCanvas("", "", 800, 600);
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("#lambda");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., ymax);
  template_h -> SetLineWidth(0);
  template_h -> Draw();

  int set_Nptx = 1000.;
  vav_kappa_0p005 -> SetNpx(set_Nptx);
  vav_kappa_0p002 -> SetNpx(set_Nptx);
  vav_kappa_0p001 -> SetNpx(set_Nptx);
  landau_c_1 -> SetNpx(set_Nptx);
  landau_c_piover2 -> SetNpx(set_Nptx);

  vav_kappa_0p005 -> SetLineColor(kBlue);
  vav_kappa_0p002 -> SetLineColor(kCyan);
  vav_kappa_0p001 -> SetLineColor(kGreen);
  vav_kappa_0p005 -> SetLineWidth(4);
  vav_kappa_0p002 -> SetLineWidth(3);
  vav_kappa_0p001 -> SetLineWidth(2);
  vav_kappa_0p002 -> SetLineStyle(7);
  vav_kappa_0p001 -> SetLineStyle(4);

  landau_c_1 -> SetLineColor(kRed);
  landau_c_piover2 -> SetLineColor(kRed);
  landau_c_piover2 -> SetLineStyle(7);

  vav_kappa_0p005 -> Draw("lsame");
  vav_kappa_0p002 -> Draw("lsame");
  vav_kappa_0p001 -> Draw("lsame");

  landau_c_1 -> Draw("lsame");
  landau_c_piover2 -> Draw("lsame");

  TLegend *l = new TLegend(0.5, 0.5, 0.85, 0.85);
  l -> SetBorderSize(0);
  l -> AddEntry(vav_kappa_0p005, "Vavilov #kappa = 0.005", "l");
  l -> AddEntry(vav_kappa_0p002, "Vavilov #kappa = 0.002", "l");
  l -> AddEntry(vav_kappa_0p001, "Vavilov #kappa = 0.001", "l");
  l -> AddEntry(landau_c_1, "Landau c = 1", "l");
  l -> AddEntry(landau_c_piover2, "Landau c = #pi/2", "l");
  l -> Draw("same");

  c -> SaveAs("./Vavilov_Landau_comparison.pdf");
  c -> Close();
}

void landau_mean_distribution(){

  TF1 *landau_c_1 = new TF1("landau_c_1", "TMath::Landau(x, 0., 1.)", -20., 200.);
  double lambda_low = -2.;
  double lambda_high = 400.;
  double lambda_scan_step = 0.01;
  int N_lambda_scan_steps = (lambda_high - lambda_low) / lambda_scan_step;

  vector<double> lambda_vec;
  vector<double> landau_mean_vec;

  for(int i = 0; i < N_lambda_scan_steps; i++){
    double this_lambda = lambda_low + lambda_scan_step * (i + 0);
    double this_mean = landau_c_1 -> Mean(-20., this_lambda);
    lambda_vec.push_back(this_lambda);
    landau_mean_vec.push_back(this_mean);
  }

  TGraph *graph = new TGraph(lambda_vec.size(), lambda_vec.data(), landau_mean_vec.data());

  double ymin = -3.;
  double ymax = 6.5;
  TCanvas *c = new TCanvas("", "", 800, 600);
  TH1D *template_h = new TH1D("", "", 1, lambda_low - 20., lambda_high + 20.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("#lambda");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("#lambda_{mean}");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(ymin, ymax * 1.2);
  template_h -> SetLineWidth(0);
  template_h -> Draw();

  graph -> Draw("psame");

  c -> SaveAs("./lambda_mean_distribution.pdf");
  
}

void lambda_max_formula_comparison(){

  double ymax = 5000.;
  
  TF1 *landau_c_1 = new TF1("landau_c_1", "TMath::Landau(x, 0., 1.)", -20., ymax);

  double lambda_low = -3.;
  double lambda_high = 8.;
  double lambda_scan_step = 0.001;
  int N_lambda_scan_steps = (lambda_high - lambda_low) / lambda_scan_step;
  
  double lambda_max_scan_low = -3.;
  double lambda_max_scan_high = ymax;
  double lambda_max_scan_step_small = 0.0001;
  double lambda_max_scan_step_large = 0.01;

  //int N_lambda_max_scan_steps = (lambda_max_scan_high - lambda_max_scan_low) / 

  vector<double> lambda_vec;
  vector<double> lambda_max_vec;

  double previoius_lambda_max = lambda_max_scan_low;
  for(int i = 0; i < N_lambda_scan_steps; i++){
    if(i % 10 == 0) cout << i << " / " << N_lambda_scan_steps << endl;
    double this_lambda = lambda_low + lambda_scan_step * (i + 0);

    double this_lambda_max = previoius_lambda_max;
    double previous_diff = 9999.;
    bool diff_increasing = false;
    while(!diff_increasing){

      double this_mean = landau_c_1 -> Mean(-20., this_lambda_max);
      double this_diff = fabs(this_lambda - this_mean);
      diff_increasing = this_diff > previous_diff;
      if(diff_increasing && previous_diff > 0.001) this_diff = false; 
      
      if(diff_increasing) cout << "this_lambda : " << this_lambda << ", this_lambda_max : " << this_lambda_max << ", this_mean : " << this_mean << ", previous_diff : " << previous_diff << ", this_diff : " << this_diff << ", diff_increasing : " << diff_increasing << endl;

      if(this_lambda_max > lambda_max_scan_high){
	cout << "Saturated, break for this_lambda : " << this_lambda << endl;
	break;
      }

      if(this_lambda_max < 30.) this_lambda_max = this_lambda_max + lambda_max_scan_step_small;
      else this_lambda_max = this_lambda_max + lambda_max_scan_step_large;
      previous_diff = this_diff;
    }

    cout << "this_lambda : " << this_lambda << ", this_lambda_max : " << this_lambda_max << endl;
    lambda_vec.push_back(this_lambda);
    lambda_max_vec.push_back(this_lambda_max);
    previoius_lambda_max = this_lambda_max;
  }

  TGraph *graph = new TGraph(lambda_vec.size(), lambda_vec.data(), lambda_max_vec.data());

  TCanvas *c = new TCanvas("", "", 800, 600);

  TH1D *template_h = new TH1D("", "", 1, lambda_low - 1., lambda_high + 1.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("#lambda");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("#lambda_{max}");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(-200., ymax * 1.2);
  template_h -> SetLineWidth(0);
  template_h -> Draw();  
  
  graph -> Draw("lpsame");
  
  //0.60715 + 1.1934 * lambda_bar + (0.67794 + 0.052382 * lambda_bar) * exp(0.94753 + 0.74442 * lambda_bar);
  TF1 *geant_lambda_max_f = new TF1("geant_lambda_max_f", "[0] + [1] * x + ([2] + [3] * x) * exp([4] + [5] * x)", lambda_low, lambda_high); 
  geant_lambda_max_f -> SetParameters(0.60715, 1.1934, 0.67794, 0.052382, 0.94753, 0.74442);
  geant_lambda_max_f -> SetLineColor(kRed);
  geant_lambda_max_f -> SetLineStyle(3);
  geant_lambda_max_f -> Draw("lsame");

  TF1 *fit_lambda_max_f = new TF1("fit_lambda_max_f", "[0] + [1] * x + ([2] + [3] * x) * exp([4] + [5] * x)", lambda_low, lambda_high);
  fit_lambda_max_f -> SetParameters(0.60715, 1.1934, 0.67794, 0.052382, 0.94753, 0.74442);
  fit_lambda_max_f -> SetLineColor(kGreen);
  fit_lambda_max_f -> SetLineStyle(7);

  graph -> Fit("fit_lambda_max_f", "RN");
  fit_lambda_max_f -> Draw("lsame");

  TLegend *l = new TLegend(0.2, 0.6, 0.6, 0.85);
  l -> SetBorderSize(0);
  l -> AddEntry(graph, "#lambda_{max}", "pl");
  l -> AddEntry(fit_lambda_max_f, "p_{0} + p_{1}#lambda + (p_{2} + p_{3}#lambda)exp(p_{4} + p_{5}#lambda)", "");
  l -> AddEntry(geant_lambda_max_f, "Traditional fit parameters", "l");
  l -> AddEntry(fit_lambda_max_f, "New fit parameters", "l");
  l -> Draw("same");

  c -> SaveAs("./lambda_max_formula_comparison.pdf");
  c -> Close();

  graph -> SetName("labmda_max_gr");
  graph -> Write();
}

void comparison(){

  out_rootfile = new TFile("./output.root", "RECREATE");
  out_rootfile -> cd();

  Landau_Vavilov_comparison();
  landau_mean_distribution();
  lambda_max_formula_comparison();

  out_rootfile -> Close();
}
	       
