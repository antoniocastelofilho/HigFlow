#include <torch/script.h> // LibTorch
#include <iostream>
#include <vector>
#include <cstring> // Para memcpy

// Variável global para manter o modelo na memória durante a simulação
torch::jit::script::Module mlp_module;
bool model_loaded = false;

// O extern "C" é mágico: ele desliga o "C++" na assinatura da função,
// permitindo que o código em C puro consiga enxergá-la e chamá-la.
extern "C" {

    // Função 1: Carrega o modelo (Chame apenas 1 vez no início do HigFlow)
    int init_mlp_model(const char* model_path) {
        try {
            mlp_module = torch::jit::load(model_path);
            mlp_module.eval(); // Garante que o Dropout/BatchNorm estão desligados
            model_loaded = true;
            return 0; // 0 significa Sucesso
        }
        catch (const c10::Error& e) {
            std::cerr << "[C++ Erro] Falha ao carregar o modelo:\n" << e.what() << std::endl;
            return -1; // -1 significa Erro
        }
    }

    // Função 2: Executa a rede (Será chamada para cada estêncil/batch)
    // Recebe o array de features do C, e um array vazio onde os pesos serão salvos
    int run_mlp_inference(const float* input_features, int max_n, int features_por_no, float* output_weights) {
        if (!model_loaded) {
            std::cerr << "[C++ Erro] O modelo não foi carregado. Chame init_mlp_model primeiro!" << std::endl;
            return -1;
        }

        int tamanho_entrada = max_n * features_por_no;

        try {
            // 1. Converte o Array do C em um Tensor do PyTorch (sem duplicar memória)
            // O formato {1, tamanho_entrada} significa: Batch Size = 1, Features = tamanho_entrada
            torch::Tensor input_tensor = torch::from_blob((void*)input_features, {1, tamanho_entrada}, torch::kFloat32);

            // 2. Prepara a entrada
            std::vector<torch::jit::IValue> inputs;
            inputs.push_back(input_tensor);

            // 3. Executa o modelo (Forward Pass)
            torch::Tensor output_tensor = mlp_module.forward(inputs).toTensor();

            // 4. Copia os resultados de volta para o array do C
            output_tensor = output_tensor.contiguous();
            std::memcpy(output_weights, output_tensor.data_ptr<float>(), max_n * sizeof(float));

            return 0; // Sucesso
        }
        catch (const c10::Error& e) {
            std::cerr << "[C++ Erro] Falha durante a inferencia:\n" << e.what() << std::endl;
            return -1;
        }
    }

} // Fim do extern "C"