#define DRAW

void solve_function()
{
  ifstream ffi;
  ffi.open("../data/Ca_Fe_modified.dat");

  vector<double> v_x;
  vector<double> v_y_Ca;
  vector<double> v_y_Fe;
  
  double a, b, c;
  while(1){ 
    ffi >> a >> b >> c;
    if(!ffi.good()) break;
    v_x.push_back(a);
    v_y_Ca.push_back(b);
    v_y_Fe.push_back(c);
  }
  ffi.close();
  
  int size = v_x.size();
  TGraph *gr_Ca = new TGraph(size);
  TGraph *gr_Fe = new TGraph(size);
  for(int i=0;i<size;i++){
    gr_Ca->SetPoint(i, v_x[i], v_y_Ca[i]);
    gr_Fe->SetPoint(i, v_x[i], v_y_Fe[i]);
  }
  
#ifdef DRAW
  TCanvas *cc = new TCanvas("cc", "cc", 0, 0, 480, 360); 
  cc->cd();
  gr_Ca->SetLineColor(2);
  gr_Fe->SetLineColor(4);
  gr_Ca->Draw("");
  gr_Fe->Draw("SAME");
#endif

  cout << "x = 1, " << gr_Ca->Eval(1.) << endl;
  cout << "x = 1, " << gr_Fe->Eval(1.) << endl;
}
