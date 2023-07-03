# Define a helper function in C++ using Rcpp
library(Rcpp)
thresholds <- seq(1000, 5000, 1000)
seeds <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
vars <- c("BMI", "AGE", "SEX", "MHABNWBC")

cppFunction('

  RObject callWithOne(Function f) {
    return f(1);
  }

')

cppFunction('
List collect_bestmodels_cpp(NumericVector thresholds, NumericVector seeds, CharacterVector vars) {
  List cDEG_bestmodel_list;
  List noncDEG_bestmodel_list;
  
  for (int i = 0; i < thresholds.size(); ++i) {
    double threshold_i = thresholds[i];
    
    for (int j = 0; j < seeds.size(); ++j) {
      double seed_i = seeds[j];
      std::string file_path = "../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti" + std::to_string(seed_i) + "FDR1.RData";
      Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
      Rcpp::Function load = base["load"];
      load(file_path);
      
      for (int k = 0; k < vars.size(); ++k) {
        std::string var_i = Rcpp::as<std::string>(vars[k]);
        List results = collect_bestmodel(threshold_i, seed_i, var_i);
        RObject callWithOne(Function collect_bestmodel) {return collect_bestmodel(threshold_i, seed_i, var_i);}
        
        
        List cDEG_bestmodel_i = results["cDEG.model.result"];
        List noncDEG_bestmodel_i = results["noncDEG.model.result"];
        
        std::string name_i = var_i + "+" + std::to_string(threshold_i) + "+" + std::to_string(seed_i);
        cDEG_bestmodel_list[name_i] = cDEG_bestmodel_i;
        noncDEG_bestmodel_list[name_i] = noncDEG_bestmodel_i;
      }
    }
  }
  
  return List::create(cDEG_bestmodel_list, noncDEG_bestmodel_list);
}')

