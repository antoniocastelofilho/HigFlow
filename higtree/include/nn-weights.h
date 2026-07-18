#ifndef NN_WEIGHTS_H
#define NN_WEIGHTS_H

#pragma push_macro("DIM")
#undef DIM
#ifdef __cplusplus
#include <torch/script.h>
#endif
#pragma pop_macro("DIM")

#include "types.h"
#include "coord.h"
#include "wls.h"

#ifdef __cplusplus
extern "C" {
#endif

   /** 
    * @brief Must executed in the beginning ofthe code.
    * Loads the mlp model to do the inference of the interpolation weights
    * @param model_path path to the .pt containing the model
   */
   int init_mlp_model(const char* model_path);

      /**
    * @brief Makes the inference
    * @param max_n Max size of an interpolation scheme
    * @param features_por_no number of feature per point (using 3 currently)
    * @param input_features a vector containing
    * @param output_weights vector containing the returned weights 
    */
   void nn_inference(int numpts, wls_item items[], int max_n, real w[]);

#ifdef __cplusplus
}
#endif
#endif // NN_WEIGHTS_H