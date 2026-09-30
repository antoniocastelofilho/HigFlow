import os, re

# Padrão para capturar os includes locais: #include "arquivo.h"
include_regex = re.compile(r'^\s*#\s*include\s*["\']([^"\']+)["\']')

print("digraph G {")
print('  node [shape=box, style=filled, color="#E8F0FE", fontname="Helvetica"];')
print('  edge [color="#4285F4"];')

for root, _, files in os.walk('./../src'):
    # Ignora pastas de build, Git ou a própria documentação do Doxygen
    if any(p in root for p in ['/build', '/.git', '/html', '/latex']):
        continue
    for file in files:
        if file.endswith(('.c', '.h', '.cpp', '.hpp')):
            file_path = os.path.join(root, file).replace('./', '')
            base_file = os.path.basename(file_path)
            
            try:
                with open(os.path.join(root, file), 'r', encoding='utf-8', errors='ignore') as f:
                    for line in f:
                        match = include_regex.match(line)
                        if match:
                            included_file = os.path.basename(match.group(1))
                            # Garante que as strings fiquem entre aspas para o Graphviz não quebrar com pontos
                            print(f'  "{base_file}" -> "{included_file}";')
            except Exception:
                pass

print("}")
