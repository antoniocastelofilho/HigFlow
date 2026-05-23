#ifndef MLP_H
#define MLP_H

/* Essa macro mágica avisa: Se quem estiver lendo isso for um compilador C++, 
   use 'extern "C"'. Se for o compilador C, ignore e leia só as funções. */
#ifdef __cplusplus
extern "C" {
#endif

// As assinaturas das nossas funções
int init_mlp_model(const char* model_path);
int run_mlp_inference(const float* input_features, int max_n, int features_por_no, float* output_weights);

#ifdef __cplusplus
}
#endif

#endif // MLP_H