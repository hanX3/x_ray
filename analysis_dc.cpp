#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>

#define LABLECALISTART "<<CALIBRATION>>"
#define LABLECALIEND "<<DATA>>"
#define LABELDATASTART "<<DATA>>"
#define LABELDATAEND "<<END>>"

#define EFFENERGY1 1.85
#define EFFENERGY2 15.

//cali data
#define ENERGYCALIPAR0 -0.0257785
#define ENERGYCALIPAR1 0.0520857
#define ENERGYCALIPAR2 3.53354e-007

//rejection info
#define AGXRAYCHANNEL 59
#define AGXRAYCOUNT 1200
#define AGXRAYSIGMA 1.88
#define AGXRAYREJECTIONCHANNELS 10

//atten eff modified info
#define ATTENEFFENERGYMIN 1.5
#define ATTENEFFENERGYMAX 15.

#define DEBUG
//#define DRAW

using namespace::std;

void ReadXRayData(map<string, map<double, double>> &m_p);
void PreAnaXRayData_1(map<string, map<double, double>> &m_p);
void PreAnaXRayData_2(map<string, map<double, double>> &m_p);
void ReadEffData(double par1[10], double par2[10], double par3[10]);
void ReadMcaData(char filename[1024], vector<double> &x, vector<double> &y);
void RejectBackground(vector<double> &x, vector<double> &y, vector<double> &yy);
void ModifiedAttenEff(double par[3], vector<double> &y, vector<double> &yy);
void GetPeakInfo(double par[3], vector<double> &v_x, vector<double> &v_y, map<double, double> &m_p);
void CaliDetectorEff(double par1[10], double par2[10], double par3[10], map<double, double> &m_p, map<double, double> &m_pp);
void CaliAttenEff(map<double, double> &m_p, map<double, double> &m_pp);
void GetElementShow(vector<string> &v_e, map<string, map<double, double>> &m_p, map<double, double> &m_pp, map<string, map<double, double>> &m_ppp, map<string, map<double, double>> &m_q, map<string, map<double, double>> &m_r);
void GetElementPercent(map<string, map<double, double>> &m_p, map<string, map<double, double>> &m_pp, map<string, map<double, double>> m_ppp, map<string, double> &m_q);
void ModifiedCaFe(map<string, double> &m_p, map<string, double> &m_q);

//
void analysis_dc()
{
  gROOT->SetBatch(0);

  map<string, map<double, double>> map_xray_data;
  ReadXRayData(map_xray_data);
  PreAnaXRayData_1(map_xray_data);
  map<string, map<double, double>>::iterator it1;
  map<double, double>::iterator it2;

  //need to choose here
  //1
  vector<string> v_element = {"20Ca", "26Fe", "82Pb", "25Mn", "24Cr", "29Cu", "30Zn", "22Ti", "14Si"};
  //2
  //v_element.clear();
  //for(it1=map_xray_data.begin();it1!=map_xray_data.end();++it1) v_element.push_back(it1->first);

  double par1_eff[10];
  double par2_eff[10];
  double par3_eff[10];
  ReadEffData(par1_eff, par2_eff, par3_eff);

  double par_cali[3];
  par_cali[0] = ENERGYCALIPAR0;
  par_cali[1] = ENERGYCALIPAR1;
  par_cali[2] = ENERGYCALIPAR2;
#ifdef DEBUG
  cout << "cali par: " << par_cali[0] << " " << par_cali[1] << " " << par_cali[2] << endl;
#endif

  char file_name_1[1024];
  char file_name_2[1024];
  sprintf(file_name_1, "../spectrum/0.4Fe3O4+0.1CaCO3-15KV-0.02mA-2.mca");
  sprintf(file_name_2, "../spectrum/0.4Fe3O4+0.1CaCO3-15KV-0.02mA-2.mca");

  //file 1 ana
  vector<double> v_x_1;
  vector<double> v_y_1;
  vector<double> v_y_rb_1;//rejection background
  vector<double> v_y_rb_ae_1;//rejection background & atten eff
  map<double, double> m_peak_1;
  map<double, double> m_peak_dector_eff_1;
  map<string, map<double, double>> m_element_show_1;
  map<string, map<double, double>> m_element_show_mca_toi_associate_1;
  map<string, map<double, double>> m_element_show_toi_mca_associate_1;

  ReadMcaData(file_name_1, v_x_1, v_y_1);
  RejectBackground(v_x_1, v_y_1, v_y_rb_1);
  ModifiedAttenEff(par_cali, v_y_rb_1, v_y_rb_ae_1);
  GetPeakInfo(par_cali, v_x_1, v_y_rb_ae_1, m_peak_1);
  CaliDetectorEff(par1_eff, par2_eff, par3_eff, m_peak_1, m_peak_dector_eff_1);

#ifdef DEBUG
  for(int i=0;i<3;i++){
    cout << "cali par " << i << " = " << par_cali[i] << endl;
  }
  for(it2=m_peak_1.begin();it2!=m_peak_1.end();++it2){
    cout << ",,, raw " << it2->first << " ==>  " << it2->second << endl; 
  }
  for(it2=m_peak_dector_eff_1.begin();it2!=m_peak_dector_eff_1.end();++it2){
    cout << ",,, after detector eff " << it2->first << " ==>  " << it2->second << endl; 
  }
#endif
  
  GetElementShow(v_element, map_xray_data, m_peak_dector_eff_1, m_element_show_1, m_element_show_mca_toi_associate_1, m_element_show_toi_mca_associate_1);
#ifdef DEBUG
  for(it1=m_element_show_1.begin();it1!=m_element_show_1.end();++it1){
    cout << "!!!... " << it1->first << " ";
    for(it2=it1->second.begin();it2!=it1->second.end();++it2){
      cout << it2->first << " " << it2->second;
    }
    cout << endl;
  }
#endif

  //file 2 ana
  vector<double> v_x_2;
  vector<double> v_y_2;
  vector<double> v_y_rb_2;//rejection background
  vector<double> v_y_rb_ae_2;//rejection background & atten eff
  map<double, double> m_peak_2;
  map<double, double> m_peak_dector_eff_2;
  map<string, map<double, double>> m_element_show_2;
  map<string, map<double, double>> m_element_show_mca_toi_associate_2;
  map<string, map<double, double>> m_element_show_toi_mca_associate_2;

  ReadMcaData(file_name_2, v_x_2, v_y_2);
  RejectBackground(v_x_2, v_y_2, v_y_rb_2);
  ModifiedAttenEff(par_cali, v_y_rb_2, v_y_rb_ae_2);
  GetPeakInfo(par_cali, v_x_2, v_y_rb_ae_2, m_peak_2);
  CaliDetectorEff(par1_eff, par2_eff, par3_eff, m_peak_2, m_peak_dector_eff_2);

#ifdef DEBUG
  for(int i=0;i<3;i++){
    cout << "cali par " << i << " = " << par_cali[i] << endl;
  }
  for(it2=m_peak_2.begin();it2!=m_peak_2.end();++it2){
    cout << ",,, raw " << it2->first << " ==>  " << it2->second << endl; 
  }
  for(it2=m_peak_dector_eff_2.begin();it2!=m_peak_dector_eff_2.end();++it2){
    cout << ",,, after detector eff " << it2->first << " ==>  " << it2->second << endl; 
  }
#endif
  
  GetElementShow(v_element, map_xray_data, m_peak_dector_eff_2, m_element_show_2, m_element_show_mca_toi_associate_2, m_element_show_toi_mca_associate_2);
#ifdef DEBUG
  for(it1=m_element_show_2.begin();it1!=m_element_show_2.end();++it1){
    cout << "!!!... " << it1->first << " ";
    for(it2=it1->second.begin();it2!=it1->second.end();++it2){
      cout << it2->first << " " << it2->second;
    }
    cout << endl;
  }
#endif

  //double check
  map<string, map<double, double>> m_element_show(m_element_show_1);
  map<string, map<double, double>> m_element_show_mca_toi_associate(m_element_show_mca_toi_associate_1);
  map<string, map<double, double>> m_element_show_toi_mca_associate(m_element_show_toi_mca_associate_1);
  for(it1=m_element_show_1.begin();it1!=m_element_show_1.end();++it1){
    if(m_element_show_2.find(it1->first) == m_element_show_2.end()){
      m_element_show.erase(it1->first);
      m_element_show_mca_toi_associate.erase(it1->first);
      m_element_show_toi_mca_associate.erase(it1->first);
    }
  }

  //
  map<string, double> m_percent_result;
  GetElementPercent(map_xray_data, m_element_show_1, m_element_show_mca_toi_associate_1, m_percent_result);

  for(map<string, double>::iterator ittt=m_percent_result.begin();ittt!=m_percent_result.end();++ittt){
    cout << "!!!... " << ittt->first << " " << ittt->second*100 << "%" << endl;
  }

  //
  map<string, double> m_percent_result_Ca_Fe_modified;
  ModifiedCaFe(m_percent_result, m_percent_result_Ca_Fe_modified);
  for(map<string, double>::iterator ittt=m_percent_result_Ca_Fe_modified.begin();ittt!=m_percent_result_Ca_Fe_modified.end();++ittt){
    cout << "!!!... " << ittt->first << " last " << ittt->second*100 << "%" << endl;
  }

}


/////
void ReadXRayData(map<string, map<double, double>> &m_p)
{
  //m_p: toi x-ray info
  ifstream file_in("../data/X-Ray.dat");
  if(!file_in){
    std::cout << "can not open X-Ray.dat file" << std::endl;
    return;
  }

  std::string line;
  std::string element;
  double energy, intensity;
  while(1){
    getline(file_in, line);
    std::istringstream iss(line);
    if(!file_in.good()) break;
    iss >> element;
    m_p.insert(pair<string, map<double, double>>(element, map<double,double>()));
    while(iss >> energy >> intensity){
      m_p[element].insert(pair<double, double>(energy, intensity));
    }
  }

  file_in.close();

  return;
}

//
//only keep the strongest peak in the x-ray data
//m_p: toi x-ray info
void PreAnaXRayData_1(map<string, map<double, double>> &m_p)
{
  map<string, map<double, double>>::iterator it;

  map<double, double> m_pp;//x-ray intensity & x-ray energy
  for(it=m_p.begin();it!=m_p.end();++it){
    m_pp.clear();
    int size = it->second.size();
    if(size > 1){
      for(map<double, double>::iterator itt=it->second.begin();itt!=it->second.end();++itt){
        m_pp.insert(pair<double, double>(itt->second, itt->first));
      }
      (it->second).clear();
      map<double, double>::reverse_iterator rit;
      rit=m_pp.rbegin();
      (it->second).insert(pair<double, double>(rit->second, rit->first));

    }else{
      continue;
    }
  }
}

//
//merge peaks energy difference smaller than 0.3keV in the x-ray data
void PreAnaXRayData_2(map<string, map<double, double>> &m_p)
{
  //m_p: toi x-ray info
  map<string, map<double, double>>::iterator it;

  map<double, int> m_diff;//need to be sorted
  for(it=m_p.begin();it!=m_p.end();++it){
    while(1){
      int size = it->second.size();
      //cout << "size  " << size << endl;
      if(size > 1){
        //cout << it->first << endl;

        vector<double> v_energy;
        //
        for(map<double, double>::iterator itt=it->second.begin();itt!=it->second.end();++itt){
          v_energy.push_back(itt->first);
        }
        sort(v_energy.begin(), v_energy.end());
        /*
        cout << "energy ";
        for(double d:v_energy){
          cout << d << " ";
        }
        cout << endl;
        */

        map<double, int> m_diff;
        for(int i=1;i<v_energy.size();i++){
          m_diff.insert(pair<double, int>(v_energy[i]-v_energy[i-1], i));
        }

        map<double, int>::iterator itu = m_diff.begin();
        //cout << " ,, " << itu->first << " " << itu->second << endl;
        if(itu->first<0.3){
          //cout << itu->first << " ==> " << itu->second << " !! " << v_energy[itu->second] << "  " << v_energy[itu->second-1] << endl; 
          //cout << " !! " << v_energy[itu->second] << " " << it->second[v_energy[itu->second]] << "  " << v_energy[itu->second-1] << " " << it->second[v_energy[itu->second-1]] << endl;
          double ea = v_energy[itu->second];
          double ia = it->second[v_energy[itu->second]];
          double eb = v_energy[itu->second-1];
          double ib = it->second[v_energy[itu->second-1]];
          double ee = ea*(ia/(ia+ib))+eb*(ib/(ia+ib));
          double ii = ia+ib;
          it->second.erase(ea);
          it->second.erase(eb);
          it->second.insert(pair<double, double>(ee, ii));
        }
        else  break;
        /*
        for(map<double, double>::iterator itt=it->second.begin();itt!=it->second.end();++itt){
          cout << " ?? " << itt->first << "  " << itt->second << endl;
        }
        */
      }
      else  break;
    }//while
  }
}

//
void ReadEffData(double par1[10], double par2[10], double par3[10])
{
  //par1, par2, par3: curve of 3 cubric function for detector efficiency
  //read effdata
  double par_temp1, par_temp2, par_temp3;
  ifstream ifs_eff;
  ifs_eff.open("../data/eff.dat");
  if(!ifs_eff){
    cout << "cannot open eff.par file !" << endl;
    return;
  }
  int kk = 0;
  while(1){
    ifs_eff >> par_temp1 >> par_temp2 >> par_temp3;
    if(!ifs_eff.good()) break;
    par1[kk] = par_temp1;
    par2[kk] = par_temp2;
    par3[kk] = par_temp3;
    kk++;
  }
  ifs_eff.close();
  /*
  for(int i=0;i<10;i++){
    cout << par1[i] << " " << par2[i] << " " << par3[i] << endl;
  }
  */
}

//
void ReadMcaData(char filename[1024], vector<double> &x, vector<double> &y)
{
  //filename: mca file
  //x:channel info
  //y:counts info
  ifstream fi;
  fi.open(filename);
  if(!fi){
    cout << "cannot open mca file ..." << endl;
    return;
  }

  string temp;

  vector<double> cali_ch;
  vector<double> cali_e;
  double cali_ch_ = 0.;
  double cali_e_ = 0.;

  while(1){
    temp.clear();
    getline(fi, temp);
    temp.erase(temp.find_last_not_of("\n\r") + 1);
    if (!fi.good())   break;

    //read cali points
    if(temp.compare(LABLECALISTART) == 0){
      getline(fi, temp);
      temp.erase(temp.find_last_not_of("\n\r") + 1);
      while(1){
        getline(fi, temp);
        temp.erase(temp.find_last_not_of("\n\r") + 1);

        if(temp.compare(LABLECALIEND) == 0){
          break;
        }

        stringstream ss;
        ss.str("");
        ss << temp;
        ss >> cali_ch_ >> cali_e_;
        cali_ch.push_back(cali_ch_);
        cali_e.push_back(cali_e_);
      }
    }

    //read data
    if(temp.compare(LABELDATASTART) == 0){
      double xx = 1;
      while(1){
        getline(fi, temp);
        temp.erase(temp.find_last_not_of("\n\r") + 1);
        if(temp.compare(LABELDATAEND) == 0){
          break;
        }
        x.push_back(xx);
        y.push_back((double)atoi(temp.c_str()));
        xx += 1.;
      }
    }
  }

  fi.close();
}

//
void RejectBackground(vector<double> &x, vector<double> &y, vector<double> &yy)
{
  //x: channel info of mca to analysis
  //y: count info of mca to analysis
  //yy: count info of mca to analysis after bacground rejection
  vector<double> v_x_back;
  vector<double> v_y_back;
  char file_name_back[1024];
  sprintf(file_name_back, "../data/background-15KV-0.02mA.mca");
  ReadMcaData(file_name_back, v_x_back, v_y_back);

  for(int i=0;i<x.size();i++){
    yy.push_back(y[i]);
    if(i>=(AGXRAYCHANNEL-AGXRAYREJECTIONCHANNELS) && i<=(AGXRAYCHANNEL+AGXRAYREJECTIONCHANNELS)){
      yy[i] = y[i] - AGXRAYCOUNT*TMath::Gaus((Double_t)i, AGXRAYCHANNEL, AGXRAYSIGMA);
    }
    else{
      yy[i] = y[i] - v_y_back[i];
    }

    if(yy[i]<0){
      yy[i] = 0;
    }
  }

#ifdef DRAW
  int bin_number = x.size();
  TCanvas *cc = new TCanvas("cc", "cc", 0, 0, 480, 360);
  TH1D *h1 = new TH1D("h1", "h_spectrum", bin_number, x[0], x[bin_number-1]);
  for(int i=0;i<bin_number;i++){
    h1->SetBinContent(x[i], y[i]);
  }
  h1->GetXaxis()->SetTitle("Channel");
  h1->GetYaxis()->SetTitle("Count");
  cc->cd();
  h1->Draw();

  int bin_number_back = v_x_back.size();
  TCanvas *cc_back = new TCanvas("cc_back", "cc_back", 500, 0, 480, 360);
  TH1D *h2 = new TH1D("h2", "h_back", bin_number_back, v_x_back[0], v_x_back[bin_number-1]);
  for(int i=0;i<bin_number_back;i++){
    h2->SetBinContent(v_x_back[i], v_y_back[i]);
  }
  h2->GetXaxis()->SetTitle("Channel");
  h2->GetYaxis()->SetTitle("Count");
  cc_back->cd();
  h2->Draw();

  TCanvas *cc3 = new TCanvas("cc3", "cc3", 1000, 0, 480, 360);
  TH1D *h3 = new TH1D("h3", "h_rejection", bin_number, x[0], x[bin_number-1]);
  h3->Add(h1, h2, 1, -1);
  cc3->cd();
  h3->Draw();

  TCanvas *cc4 = new TCanvas("cc4", "cc4", 0, 400, 480, 360);
  TH1D *h4 = new TH1D("h4", "h_rejection", bin_number, x[0], x[bin_number-1]);
  for(int i=0;i<bin_number;i++){
    h4->SetBinContent(x[i], yy[i]);
  }
  cc4->cd();
  h4->Draw();
#endif
}

//
void ModifiedAttenEff(double par[3], vector<double> &y, vector<double> &yy)
{
  //par[3]: cali info
  //y: peak info from mca spectrum after background rejection
  //yy: peak info from mca spectrum after background rejection and atten effective modified
  ifstream fi;
  fi.open("../data/gamma_attenuation.dat");
  if(!fi){
    cout << "cannot open gamma_attenuation_data.txt file !" << endl;
    return;
  }
  vector<double> v_x_ray_e;
  vector<double> v_x_ray_eff;

  double a, b;
  while(1){
    fi >> a >> b;
    if(!fi.good()) break;
    v_x_ray_e.push_back(a);
    v_x_ray_eff.push_back(b);
  }
  fi.close();

  int size = v_x_ray_e.size();
  TGraph *gr = new TGraph(size);
  for(int i=0;i<size;i++){
    gr->SetPoint(i, v_x_ray_e[i], v_x_ray_eff[i]);
  }

  double e = 0.;
  for(int i=0;i<y.size();i++){
    yy.push_back(y[i]);
    e = par[0]+par[1]*i+par[2]*i*i;
    if(e>1.5 && e<15.){
      yy[i] = y[i]/gr->Eval(e);
    }
  }

#ifdef DRAW
  int bin_number = y.size();
  TCanvas *cc_mae = new TCanvas("cc_mae", "cc_mae", 0, 800, 480, 360);
  TH1D *h_mae = new TH1D("h_mae", "h_rejection_atten_eff", bin_number, par[0], par[0]+par[1]*(bin_number-1)+par[2]*(bin_number-1)*(bin_number-1));
  for(int i=0;i<bin_number;i++){
    h_mae->SetBinContent(i, yy[i]);
  }
  cc_mae->cd();
  h_mae->Draw();
#endif
}

//
void GetPeakInfo(double par[3], vector<double> &x, vector<double> &y, map<double, double> &m_p)
{
  //par: calibration par
  //x: channel info
  //y: counts info
  //m_p: peak info from mca spectrum
  //histogram
  int bin_number = x.size();
  TH1D* h0_gpi = new TH1D("h0_gpi", "raw data", bin_number, x[0], x[bin_number-1]);

  for(int i=0;i<bin_number;i++){
    h0_gpi->SetBinContent(x[i], y[i]);
  }
  h0_gpi->GetXaxis()->SetTitle("Channel");
  h0_gpi->GetYaxis()->SetTitle("Count");

  //find peak
  TSpectrum *spec = new TSpectrum(20);
  int nfound = spec->Search(h0_gpi, 1, "", 0.01);
  TH1 *hb = spec->Background(h0_gpi, 20, "Compton BackSmoothing13 same");

  TF1 *tf_gaus = new TF1("tf_gaus", "gaus");

  double *xpeaks;
  xpeaks = spec->GetPositionX();
  for(int i=0;i<nfound;i++){
    double xp = xpeaks[i];
    if((par[0]+par[1]*xp+par[2]*xp*xp) < 1.5)  continue;

    h0_gpi->Fit(tf_gaus, "QN0", "", 0.97*xp, 1.03*xp);
    map<double, double> m_single_peak;
    map<double, double>::iterator it;
    for(int j=0.9*xp;j<1.1*xp;j++){
      m_single_peak.insert(pair<double, double>(h0_gpi->GetBinCenter(j), h0_gpi->GetBinContent(j)));
    }

    map<double, double> m_find_max_x;
    for(int j=xp-2;j<=xp+2;j++){
      m_find_max_x.insert(pair<double, double>(h0_gpi->GetBinContent(j), h0_gpi->GetBinCenter(j)));
    }

    map<double, double>::reverse_iterator rit = m_find_max_x.rbegin();

    double x_left_stop = 0;
    double x_right_stop = 0;
    //move left
    it = m_single_peak.find(rit->second);
    double y_max = it->second;
    while(1){
       if((it->first-(int)(0.9*xp) < 1.)){
         x_left_stop = it->first;
         break;
       }
      --it;
      if(it->second < y_max){
        y_max = it->second;
      }
      else{
        ++it;
        x_left_stop = it->first;
        break;
      }
    }
    //move right
    it = m_single_peak.find(rit->second);
    y_max = it->second;
    while(1){
       if(((int)(1.1*xp)-it->first < 1.)){
         x_right_stop = it->first;
         break;
       }
      ++it;
      if(it->second < y_max){
        y_max = it->second;
      }
      else{
        --it;
        x_right_stop = it->first;
        break;
      }
    }
    double area = h0_gpi->Integral(x_left_stop, x_right_stop) - hb->Integral(x_left_stop, x_right_stop);
    m_p.insert(pair<double, double>(par[0]+par[1]*xp+par[2]*xp*xp, area));
#ifdef DEBUG
    cout << "x_left_stop = " << x_left_stop << " x_right_stop " << x_right_stop << endl;
    cout << "area = " << area << endl;
#endif
  }

#ifdef DRAW
  TCanvas *cc_gpi = new TCanvas("cc_gpi", "cc_gpi", 0, 0, 480, 360);
  cc_gpi->cd();
  h0_gpi->Draw();
#endif
}

//
void CaliDetectorEff(double par1[10], double par2[10], double par3[10], map<double, double> &m_p, map<double, double> &m_pp)
{
  //par1, par2, par3: detector eff par
  //m_p: peak info from mca spectrum
  //m_pp: peak info from mca spectrum after detector efficiency calibration
  double eff = 0.;
  map<double, double>::iterator it;
  for(it=m_p.begin();it!=m_p.end();++it){
    eff = 0.;
    for(int i=0;i<10;i++){
      if(it->first<EFFENERGY1){
        eff += par1[i]*pow(it->first, (double)i);
      }
      else if(it->first>EFFENERGY2){
        eff += par3[i]*pow(it->first, (double)i);
      }else{
        eff += par2[i]*pow(it->first, (double)i);
      }
    }
    m_pp.insert(pair<double, double>(it->first, it->second/eff));
  }
}

//
void CaliAttenEff(map<double, double> &m_p, map<double, double> &m_pp)
{
  //m_p: peak info from mca spectrum after detector efficiency calibration
  //m_pp: peak info from mca spectrum after detector&atten efficiency calibration
  ifstream fi;
  fi.open("../data/gamma_attenuation.dat");
  if(!fi){
    cout << "cannot open gamma_attenuation_data.txt file !" << endl;
    return;
  }
  vector<double> v_x_ray_e;
  vector<double> v_x_ray_eff;

  double a, b;
  while(1){
    fi >> a >> b;
    if(!fi.good()) break;
    v_x_ray_e.push_back(a);
    v_x_ray_eff.push_back(b);
  }
  fi.close();

  int size = v_x_ray_e.size();
  TGraph *gr = new TGraph(size);
  for(int i=0;i<size;i++){
    gr->SetPoint(i, v_x_ray_e[i], v_x_ray_eff[i]);
  }

  map<double, double>::iterator it;
  for(it=m_p.begin();it!=m_p.end();++it){
    m_pp.insert(pair<double, double>(it->first, it->second/gr->Eval(it->first)));
  }
}

//
void GetElementShow(vector<string> &v_e, map<string, map<double, double>> &m_p, map<double, double> &m_pp, map<string, map<double, double>> &m_ppp, map<string, map<double, double>> &m_q, map<string, map<double, double>> &m_r)
{
  //v_e: the elements wants to caliculate
  //m_p: x-ray data from toi
  //m_pp: peak info from mca spectrum after detector&atten efficiency calibration
  //m_ppp: elements show in the mca spectrum, like <20Ca, <peak energy, peak area>, peak energy&area both from mca spectrum
  //m_q: associate peak energy from peak from mca spectrum & toi
  //m_r: associate peak energy from toi & peak from mca spectrum
  map<string, map<double, double>>::iterator itt1;
  map<double, double>::iterator itt2;

  int size = v_e.size();
  for(int i=0;i<size;i++){
    itt1 = m_p.find(v_e[i]);
#ifdef DEBUG
    cout << "// analysis " << v_e[i] << " element." << endl;
#endif
    m_ppp.insert(pair<string, map<double, double>>(v_e[i], map<double, double>()));
    m_q.insert(pair<string, map<double, double>>(v_e[i], map<double, double>()));
    m_r.insert(pair<string, map<double, double>>(v_e[i], map<double, double>()));
    
#ifdef DEBUG
    for(map<double, double>::iterator ittt=m_pp.begin();ittt!=m_pp.end();++ittt){
      cout << "//// " << ittt->first << endl;
    }
#endif
    for(itt2=itt1->second.begin();itt2!=itt1->second.end();++itt2){
#ifdef DEBUG
      cout << "energy /// " << itt2->first << endl;
#endif
      map<double, double>::iterator ittt=m_pp.begin();
      double a = ittt->first;
      double d = abs(a-itt2->first);
      while(1){
        ++ittt;
        a = ittt->first;
        if(d < abs(a-itt2->first)){
          --ittt;
          break;
        }
        else{
          d = abs(a-itt2->first);
        }
      }
      //check real distance
      if(abs(d)>0.2){
#ifdef DEBUG
        cout << "!! can not find " << v_e[i] << " " << itt2->first << " keV X-Ray" << endl;
#endif
        continue;
      }
      //save the peak
      m_ppp[v_e[i]].insert(pair<double, double>(ittt->first, ittt->second));
      m_q[v_e[i]].insert(pair<double, double>(ittt->first, itt2->first));
      m_r[v_e[i]].insert(pair<double, double>(itt2->first, ittt->first));
    }
    if(m_ppp[v_e[i]].size()==0){
      m_ppp.erase(v_e[i]);
    }
    if(m_q[v_e[i]].size()==0){
      m_q.erase(v_e[i]);
    }
    if(m_r[v_e[i]].size()==0){
      m_r.erase(v_e[i]);
    }
  }

#ifdef DEBUG
  //cout result
  for(map<string, map<double, double>>::iterator itttt=m_ppp.begin();itttt!=m_ppp.end();++itttt){
    cout << "before check<<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  " << itt2->first << "  " << itt2->second;
    }
    cout << endl;
  }
  for(map<string, map<double, double>>::iterator itttt=m_q.begin();itttt!=m_q.end();++itttt){
    cout << "before <<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  mca " << itt2->first << "keV &  toi " << itt2->second;
    }
    cout << endl;
  }
  for(map<string, map<double, double>>::iterator itttt=m_r.begin();itttt!=m_r.end();++itttt){
    cout << "before <<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  toi " << itt2->first << "keV &  mca " << itt2->second;
    }
    cout << endl;
  }
#endif

  map<double, string> mm_p;
  for(itt2=m_pp.begin();itt2!=m_pp.end();++itt2){
#ifdef DEBUG
    cout << "///......  ana  " << itt2->first << endl;
#endif
    mm_p.clear();
    for(map<string, map<double, double>>::iterator ittt=m_ppp.begin();ittt!=m_ppp.end();++ittt){
      for(map<double, double>::iterator itttt=ittt->second.begin();itttt!=ittt->second.end();++itttt){
        if(itt2->first == itttt->first){
#ifdef DEBUG
          cout << "///...  " << ittt->first << " " << itt2->first << "  " << itttt->first << "  " << m_q[ittt->first][itt2->first] << endl;
#endif
          mm_p.insert(pair<double, string>(abs(m_q[ittt->first][itt2->first]-itttt->first), ittt->first));
        }
      }
    }
#ifdef DEBUG
    for(map<double, string>::iterator ittt=mm_p.begin();ittt!=mm_p.end();++ittt){
      cout << "///...  " << ittt->first << "  " << ittt->second << endl;
    }
#endif
    if(mm_p.size()<2) continue;
    map<double, string>::iterator ittt=mm_p.begin();
    while(1){
      ++ittt;
      if(ittt==mm_p.end()) break;
      m_ppp.erase(ittt->second);
      m_q.erase(ittt->second);
      m_r.erase(ittt->second);
    }
  }

#ifdef DEBUG
  for(map<string, map<double, double>>::iterator itttt=m_ppp.begin();itttt!=m_ppp.end();++itttt){
    cout << "after check<<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  " << itt2->first << "  " << itt2->second;
    }
    cout << endl;
  }
  for(map<string, map<double, double>>::iterator itttt=m_q.begin();itttt!=m_q.end();++itttt){
    cout << "after <<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  mca " << itt2->first << "keV &  toi " << itt2->second;
    }
    cout << endl;
  }
  for(map<string, map<double, double>>::iterator itttt=m_r.begin();itttt!=m_r.end();++itttt){
    cout << "after <<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  toi " << itt2->first << "keV &  mca " << itt2->second;
    }
    cout << endl;
  }
#endif
}

//
void GetElementPercent(map<string, map<double, double>> &m_p, map<string, map<double, double>> &m_pp, map<string, map<double, double>> m_ppp, map<string, double> &m_q)
{
  //m_p: x-ray data from toi
  //m_pp: elements show in the mca spectrum
  //m_ppp: associate peak energy from peak from mca spectrum & toi 

  map<string, map<double, double>>::iterator it;
  double sum = 0.;
  for(it=m_pp.begin();it!=m_pp.end();++it){
    map<double, double>::iterator itt = it->second.begin();
    sum += (itt->second)/(m_p[it->first][m_ppp[it->first][itt->first]])*100;
  }

  for(it=m_pp.begin();it!=m_pp.end();++it){
    map<double, double>::iterator itt = it->second.begin();
    m_q[it->first] = (itt->second)/(m_p[it->first][m_ppp[it->first][itt->first]])*100/sum;
  }

}

//
void ModifiedCaFe(map<string, double> &m_p, map<string, double> &m_q)
{
  //m_p: element percent result
  //m_q: element percent result after modified Ca_Fe

  if(m_p.size()!=2){
    for(map<string, double>::iterator ittt=m_p.begin();ittt!=m_p.end();++ittt){
      m_q[ittt->first] = ittt->second;
    }
    return;
  }else{
    if((m_p.find("Ca")==m_p.end()) || (m_p.find("Fe")==m_p.end())){
      for(map<string, double>::iterator ittt=m_p.begin();ittt!=m_p.end();++ittt){
        m_q[ittt->first] = ittt->second;
      }
      return;
    }
    else{
      double k = m_p["Ca"]/m_p["Fe"];

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

      //cout << "aaa = " << v_x[0] << endl;
      //cout << "bbb = " << v_x[size-1] << endl;

      if(k>=v_x[0] && k<=v_x[size-1]){
        for(map<string, double>::iterator ittt=m_p.begin();ittt!=m_p.end();++ittt){
          if(ittt->first.compare("Ca")==0){
            m_q[ittt->first] = gr_Ca->Eval(k)*ittt->second;
          }
          else{
            m_q[ittt->first] = gr_Fe->Eval(k)*ittt->second;
          }
        }

        double a = m_q["Ca"];
        double b = m_q["Fe"];

        m_q["Ca"] = a/(a+b);
        m_q["Fe"] = b/(a+b);
      }
      else{
        for(map<string, double>::iterator ittt=m_p.begin();ittt!=m_p.end();++ittt){
          m_q[ittt->first] = ittt->second;
        }
      }

      return;
    }
  }
}
