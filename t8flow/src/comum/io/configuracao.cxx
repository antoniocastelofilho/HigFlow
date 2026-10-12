#include "configuracao.hxx"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include <libfyaml.h>

namespace t8flow {
namespace {

//! Percorre a arvore do documento e achata em `caminho.pontilhado -> texto`.
//
// SO' AS FOLHAS ESCALARES viram chave.  Sequencias viram `chave.0`, `chave.1`,
// o que preserva a ordem e permite ler uma lista sem um tipo proprio.
void achata(struct fy_node *no, const std::string &prefixo,
            const std::string &arquivo,
            std::map<std::string, std::string> *val,
            std::map<std::string, std::string> *de)
{
  if (no == nullptr) return;

  switch (fy_node_get_type(no)) {
    case FYNT_SCALAR: {
      size_t n = 0;
      const char *t = fy_node_get_scalar(no, &n);
      const std::string texto(t ? t : "", n);
      const auto it = val->find(prefixo);
      if (it != val->end() && it->second != texto)
        throw std::runtime_error(
            "Configuracao: a chave \"" + prefixo + "\" aparece duas vezes com "
            "valores diferentes -- \"" + it->second + "\" em " +
            (*de)[prefixo] + " e \"" + texto + "\" em " + arquivo +
            ".  Qual vale nao pode depender da ordem de leitura.");
      (*val)[prefixo] = texto;
      (*de)[prefixo] = arquivo;
      break;
    }
    case FYNT_MAPPING: {
      void *iter = nullptr;
      struct fy_node_pair *par = nullptr;
      while ((par = fy_node_mapping_iterate(no, &iter)) != nullptr) {
        struct fy_node *ch = fy_node_pair_key(par);
        size_t n = 0;
        const char *t = fy_node_get_scalar(ch, &n);
        if (t == nullptr) continue;          // chave nao escalar: nao ocorre no esquema
        const std::string nome(t, n);
        achata(fy_node_pair_value(par),
               prefixo.empty() ? nome : prefixo + "." + nome, arquivo, val, de);
      }
      break;
    }
    case FYNT_SEQUENCE: {
      void *iter = nullptr;
      struct fy_node *item = nullptr;
      int i = 0;
      while ((item = fy_node_sequence_iterate(no, &iter)) != nullptr)
        achata(item, prefixo + "." + std::to_string(i++), arquivo, val, de);
      break;
    }
    default:
      break;
  }
}

bool arquivo_existe(const std::string &caminho)
{
  std::FILE *f = std::fopen(caminho.c_str(), "r");
  if (f == nullptr) return false;
  std::fclose(f);
  return true;
}

}  // namespace

void Configuracao::junta(const std::string &caminho)
{
  struct fy_document *doc = fy_document_build_from_file(nullptr, caminho.c_str());
  if (doc == nullptr)
    throw std::runtime_error("Configuracao: nao consegui ler \"" + caminho +
                             "\" -- " + std::strerror(errno));
  achata(fy_document_root(doc), "", caminho, &val_, &de_);
  fy_document_destroy(doc);
}

void Configuracao::junta_caso(const std::string &prefixo)
{
  // OS QUATRO DO HIGFLOW, na ordem em que ele os nomeia.  Ausente e' legitimo
  // -- um caso sem condicao inicial comeca do repouso --, ilegivel nao e'.
  static const char *sufixos[] = {".par.contr.yaml", ".domain.yaml",
                                  ".bc.yaml", ".init.yaml"};
  int achados = 0;
  for (const char *s : sufixos) {
    const std::string caminho = prefixo + s;
    if (!arquivo_existe(caminho)) continue;
    junta(caminho);
    achados++;
  }
  if (achados == 0)
    throw std::runtime_error(
        "Configuracao: nenhum dos quatro arquivos de \"" + prefixo +
        "\" existe (.par.contr.yaml, .domain.yaml, .bc.yaml, .init.yaml).  "
        "Um caso sem nenhum deles nao e' um caso.");
}

const std::string *Configuracao::acha(const std::string &chave) const
{
  const auto it = val_.find(chave);
  if (it == val_.end()) return nullptr;
  lidas_.insert(chave);
  return &it->second;
}

void Configuracao::faltou(const std::string &chave) const
{
  // A MENSAGEM TRAZ AS CHAVES VIZINHAS, porque o erro quase sempre e' de
  // digitacao ou de nivel trocado, e ver o bloco ao lado resolve na hora.
  std::string perto;
  const size_t ponto = chave.rfind('.');
  const std::string bloco = (ponto == std::string::npos) ? "" : chave.substr(0, ponto + 1);
  int n = 0;
  for (const auto &kv : val_)
    if (!bloco.empty() && kv.first.compare(0, bloco.size(), bloco) == 0 && n < 8) {
      perto += (n++ ? ", " : "") + kv.first;
    }
  throw std::runtime_error(
      "Configuracao: falta a chave obrigatoria \"" + chave + "\"" +
      (perto.empty() ? std::string(". O bloco nao existe em nenhum arquivo lido.")
                     : ".  No mesmo bloco existem: " + perto + "."));
}

bool Configuracao::tem(const std::string &chave) const
{
  return val_.find(chave) != val_.end();
}

double Configuracao::real(const std::string &chave) const
{
  const std::string *s = acha(chave);
  if (s == nullptr) faltou(chave);
  char *fim = nullptr;
  const double v = std::strtod(s->c_str(), &fim);
  // O TIPO E' COBRADO: `dt: muito pouco` nao e' zero, e' erro.
  while (fim != nullptr && *fim == ' ') fim++;
  if (fim == s->c_str() || (fim != nullptr && *fim != '\0'))
    throw std::runtime_error("Configuracao: \"" + chave + "\" vale \"" + *s +
                             "\", que nao e' um numero real (" + de_.at(chave) + ")");
  return v;
}

long Configuracao::inteiro(const std::string &chave) const
{
  const std::string *s = acha(chave);
  if (s == nullptr) faltou(chave);
  char *fim = nullptr;
  const long v = std::strtol(s->c_str(), &fim, 10);
  while (fim != nullptr && *fim == ' ') fim++;
  if (fim == s->c_str() || (fim != nullptr && *fim != '\0'))
    throw std::runtime_error("Configuracao: \"" + chave + "\" vale \"" + *s +
                             "\", que nao e' um numero inteiro (" + de_.at(chave) + ")");
  return v;
}

bool Configuracao::booleano(const std::string &chave) const
{
  const std::string *s = acha(chave);
  if (s == nullptr) faltou(chave);
  if (*s == "true" || *s == "True" || *s == "yes" || *s == "1") return true;
  if (*s == "false" || *s == "False" || *s == "no" || *s == "0") return false;
  throw std::runtime_error("Configuracao: \"" + chave + "\" vale \"" + *s +
                           "\", que nao e' true nem false (" + de_.at(chave) + ")");
}

std::string Configuracao::texto(const std::string &chave) const
{
  const std::string *s = acha(chave);
  if (s == nullptr) faltou(chave);
  return *s;
}

double Configuracao::real(const std::string &chave, double padrao) const
{
  return tem(chave) ? real(chave) : padrao;
}

long Configuracao::inteiro(const std::string &chave, long padrao) const
{
  return tem(chave) ? inteiro(chave) : padrao;
}

bool Configuracao::booleano(const std::string &chave, bool padrao) const
{
  return tem(chave) ? booleano(chave) : padrao;
}

std::string Configuracao::texto(const std::string &chave,
                                const std::string &padrao) const
{
  return tem(chave) ? texto(chave) : padrao;
}

int Configuracao::opcao(const std::string &chave,
                        const std::vector<std::string> &validos) const
{
  const std::string v = texto(chave);
  for (size_t i = 0; i < validos.size(); i++)
    if (validos[i] == v) return static_cast<int>(i);

  // A MENSAGEM TRAZ A LISTA.  Sem ela, "valor invalido" obriga quem le' a
  // procurar no codigo quais sao os validos -- e foi para nao cair no padrao
  // em silencio que esta funcao existe.
  std::string lista;
  for (size_t i = 0; i < validos.size(); i++) lista += (i ? ", " : "") + validos[i];
  throw std::runtime_error("Configuracao: \"" + chave + "\" vale \"" + v +
                           "\", que nao e' uma opcao valida.  Validas: " + lista +
                           " (" + de_.at(chave) + ")");
}

int Configuracao::opcao(const std::string &chave,
                        const std::vector<std::string> &validos, int padrao) const
{
  return tem(chave) ? opcao(chave, validos) : padrao;
}

std::vector<std::string> Configuracao::nao_lidas() const
{
  std::vector<std::string> fora;
  for (const auto &kv : val_)
    if (lidas_.find(kv.first) == lidas_.end()) fora.push_back(kv.first);
  return fora;
}

std::vector<std::string> Configuracao::chaves() const
{
  std::vector<std::string> todas;
  todas.reserve(val_.size());
  for (const auto &kv : val_) todas.push_back(kv.first);
  return todas;
}

const std::string &Configuracao::origem(const std::string &chave) const
{
  const auto it = de_.find(chave);
  if (it == de_.end())
    throw std::out_of_range("Configuracao::origem: \"" + chave + "\" nao foi lida");
  return it->second;
}

}  // namespace t8flow
