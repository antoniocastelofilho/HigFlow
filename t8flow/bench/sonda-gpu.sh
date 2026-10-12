#!/bin/sh
# SONDA DE GPU: o que esta maquina oferece para calcular em placa.
#
# Existe porque a escolha do back-end de GPU do t8Flow depende do que cada
# maquina TEM, e isso varia -- a decisao nao pode sair de suposicao sobre o parque.
# Rode em cada maquina candidata e junte as saidas; o resultado e' o inventario.
#
# Nao instala nada, nao exige privilegio, e diz o que FALTA quando algo falta --
# "OpenCL ausente" e "OpenCL presente sem driver registrado" sao estados
# diferentes, e confundi-los custa um dia.
#
# Uso:  sh bench/sonda-gpu.sh

echo "=== maquina ==="
uname -srm
echo "nucleos logicos:  $(nproc 2>/dev/null || echo '?')"
if command -v lscpu >/dev/null 2>&1; then
    cores=$(lscpu | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2);print $2}')
    socks=$(lscpu | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2);print $2}')
    thr=$(lscpu | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2);print $2}')
    [ -n "$cores" ] && echo "nucleos fisicos:  $((cores * socks))  (${thr} thread(s) por nucleo)"
    echo "  -- o nproc conta THREADS; medida de escalabilidade usa nucleos FISICOS"
fi

echo
echo "=== placas ==="
if command -v lspci >/dev/null 2>&1; then
    lspci 2>/dev/null | grep -iE 'vga|3d controller|display' || echo "  nenhuma placa listada"
else
    echo "  lspci ausente -- tente 'ls /sys/class/drm'"
fi

echo
echo "=== NVIDIA ==="
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv 2>/dev/null \
        || echo "  nvidia-smi presente mas falhou (driver nao carregado?)"
    command -v nvcc >/dev/null 2>&1 && echo "  nvcc: $(nvcc --version 2>/dev/null | tail -1)" \
                                    || echo "  nvcc AUSENTE (CUDA toolkit nao instalado)"
else
    echo "  nvidia-smi ausente"
fi

echo
echo "=== AMD ==="
command -v rocminfo >/dev/null 2>&1 && echo "  rocminfo presente" || echo "  rocminfo ausente"
command -v hipcc   >/dev/null 2>&1 && echo "  hipcc presente"   || echo "  hipcc AUSENTE (sem compilador HIP)"
for so in /usr/lib/x86_64-linux-gnu/libamdhip64.so.* /opt/rocm*/lib/libamdhip64.so*; do
    [ -e "$so" ] && echo "  runtime HIP: $so"
done
if [ -e /dev/kfd ]; then
    echo "  /dev/kfd presente  ($(stat -c '%U:%G %A' /dev/kfd 2>/dev/null))"
    # ACESSO, nao existencia: o no' de compute e' do grupo render, e usuario fora
    # dele nao enumera a GPU -- falha que se apresenta como "nenhum dispositivo".
    if [ -r /dev/kfd ] && [ -w /dev/kfd ]; then
        echo "  acesso ao /dev/kfd: OK"
    else
        echo "  acesso ao /dev/kfd: NEGADO -- o usuario precisa do grupo"
        echo "    $(stat -c '%G' /dev/kfd 2>/dev/null) (grupos atuais: $(id -nG))"
    fi
else
    echo "  /dev/kfd AUSENTE -- sem ele o ROCm/HIP nao ve' a GPU"
fi

echo
echo "=== OpenCL ==="
loader=""
for so in /usr/lib/x86_64-linux-gnu/libOpenCL.so.1 /usr/lib64/libOpenCL.so.1; do
    [ -e "$so" ] && loader="$so"
done
if [ -n "$loader" ]; then
    echo "  loader: $loader"
else
    echo "  loader AUSENTE (instale ocl-icd-libopencl1 ou equivalente)"
fi
icds=$(ls /etc/OpenCL/vendors/*.icd 2>/dev/null)
if [ -n "$icds" ]; then
    echo "  drivers registrados:"
    for f in $icds; do echo "    $(basename "$f") -> $(cat "$f" 2>/dev/null)"; done
else
    echo "  NENHUM driver registrado em /etc/OpenCL/vendors/"
    echo "    o loader sozinho nao calcula nada: sem ICD, clGetPlatformIDs devolve 0"
fi
if command -v clinfo >/dev/null 2>&1; then
    echo "  clinfo:"
    clinfo 2>/dev/null | grep -iE 'platform name|device name|device type|opencl c version' \
        | sed 's/^/    /' | head -12
else
    echo "  clinfo ausente -- instale para enumerar de verdade"
fi

echo
echo "=== Vulkan (alternativa de compute portavel) ==="
ls /usr/share/vulkan/icd.d/*.json >/dev/null 2>&1 \
    && echo "  ICDs: $(ls /usr/share/vulkan/icd.d/ | tr '\n' ' ')" \
    || echo "  nenhum ICD Vulkan"
command -v vulkaninfo >/dev/null 2>&1 && echo "  vulkaninfo presente" || echo "  vulkaninfo ausente"

echo
echo "=== resumo para a decisao ==="
echo "O que o t8Flow precisa de uma maquina para usar GPU no lote de stencils:"
echo "  1. uma placa de compute (NVIDIA, AMD ou Intel);"
echo "  2. um runtime que o codigo possa chamar -- ICD OpenCL registrado, ou"
echo "     CUDA/HIP com compilador;"
echo "  3. ACESSO ao no' de compute (grupo render para /dev/kfd no AMD)."
echo "Faltando qualquer um dos tres, o caminho de CPU continua sendo o unico, e e'"
echo "por isso que ele e' o oraculo e nao o plano B."
