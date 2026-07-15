#include "nn-weights.h"
#include <iostream>
#include <vector>
#include <cstring>

// global  variables in pure C++
torch::jit::script::Module mlp_module;
bool model_loaded = false;

std::vector<real> buffer_global;
bool buffer_alocado = false;

torch::Tensor _data_to_tensor(int numpts, wls_item items[], int max_n, int n_features) {
    int tensor_size = max_n * n_features;
    
    // allocs memory only in the first time it runs
    if (!buffer_alocado) {
        buffer_global.resize(tensor_size, 0.0);
        buffer_alocado = true;
    } else {
        // Zera o buffer rapidamente para o zero-padding
        std::fill(buffer_global.begin(), buffer_global.end(), 0.0);
    }

    for(int i = 0; i < max_n; i++) {
        if(i < numpts) {
            for(int j = 0; j < DIM; j++) buffer_global[3*i + j] = items[i].x[j];
            buffer_global[3*i + DIM] = 1.0 / (items[i].dist * items[i].dist + EPSMACH);
        }
    }

    // creates the tensor pointing to the global buffer without cloning
    return torch::from_blob(buffer_global.data(), {1, tensor_size}, torch::kFloat64);
}

int _run_mlp_inference(torch::Tensor input_tensor, int max_n, real* output_weights) {
    if (!model_loaded) {
        std::cerr << "[C++ Erro] O modelo não foi carregado. Chame init_mlp_model primeiro!" << std::endl;
        return -1;
    }

    try {
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(input_tensor);

        torch::Tensor output_tensor = mlp_module.forward(inputs).toTensor();
        output_tensor = output_tensor.contiguous();

        std::memcpy(output_weights, output_tensor.data_ptr<real>(), max_n * sizeof(real));

        return 0;
    }
    catch (const c10::Error& e) {
        std::cerr << "[C++ Erro] Falha durante a inferencia:\n" << e.what() << std::endl;
        return -1;
    }
}


// =======================================================
//                      C INTERFACE 
// =======================================================
extern "C" {

    int init_mlp_model(const char* model_path) {
        try {
            // kCPU adicionado para evitar aquele erro de compatibilidade de GPUs
            mlp_module = torch::jit::load(model_path, torch::kCPU);
            mlp_module.eval();
            model_loaded = true;
            return 0;
        }
        catch (const c10::Error& e) {
            std::cerr << "[C++ Erro] Falha ao carregar o modelo:\n" << e.what() << std::endl;
            return -1; 
        }
    }

    void nn_inference(int numpts, wls_item items[], int max_n, real w[]) {
        // Considerando as suas features (coordenadas x,y + dist), n_features = DIM + 1
        int n_features = DIM + 1; 
        
        // 1. creates the Tensor from the HiG-Flow data
        torch::Tensor input = _data_to_tensor(numpts, items, max_n, n_features);
        
        // 2. runs the inference
        _run_mlp_inference(input, max_n, w);
    }

} // Fim do extern "C"