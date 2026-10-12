// O ARRANJO DESLOCADO DO T8FLOW 2D: onde cada propriedade MORA.
//
//   escalares  (pressao, fracao volumetrica, densidade, viscosidade)  CENTRO
//   tensoriais (tensor de tensoes, gradiente de velocidade)           CENTRO
//   vetoriais  (velocidade, fluxo)                                    FACETA
//
// NAO E' CONVENCAO ARBITRARIA.  Guardar a velocidade nas facetas e a pressao nos
// centros e' o arranjo deslocado (MAC): a divergencia de uma celula sai da soma
// das velocidades das SUAS facetas, e o gradiente de pressao numa faceta sai da
// diferenca entre os DOIS centros que ela separa -- ambos com estencil de duas
// celulas e sem interpolacao.  No arranjo colocado (tudo no centro) o gradiente
// precisa de tres celulas e admite um modo de pressao em tabuleiro que a
// divergencia nao ve.  E' o mesmo arranjo do HiGFlow, onde a velocidade vive em
// `dpu[dim]` sobre um `sim_facet_domain` POR DIRECAO e o tensor vive em
// `dpDu[i][j]`, componente por componente, nos centros.
//
// A POSSE DAS FACETAS, que e' a parte delicada e a unica decisao real aqui:
//
//   **a faceta da direcao d e' a face "menos" (face 2d) de EXATAMENTE UMA folha.**
//
// Uma faceta interna entre A (antes) e B (depois) e' a face "mais" de A e a face
// "menos" de B, logo pertence a B -- sem ambiguidade e sem negociacao.  As faces
// "mais" que sao CONTORNO DO DOMINIO nao sao face "menos" de ninguem, e por isso
// ganham faceta propria.  Em malha uniforme N x N isso da' N(N+1) facetas por
// direcao, que e' a contagem classica do MAC.
//
// E A REGRA SOBREVIVE AO REFINO 2:1 sem emenda: a face "mais" de uma celula
// grossa encostada em duas finas NAO vira uma faceta -- as duas finas possuem as
// suas faces "menos", e o lado fino fica com DOIS valores de velocidade.  E' o que
// o metodo pede ali (o fluxo pela face grossa e' a soma dos dois), e e' por isso
// que a convencao escolhida foi a face "menos" e nao a "mais": a "mais" daria uma
// faceta na celula grossa e duas nas finas para a MESMA interface.
//
// LIMITE DE HOJE: np=1.  "Face sem vizinho LOCAL" e' o criterio de contorno usado
// aqui, e sob particionamento ele passa a incluir face de FRANJA -- a distincao
// exige a camada de ghost, que o t8Flow ainda nao tem.  Esta' afirmado no teste e
// dito aqui para nao ser suposto.
#ifndef T8FLOW_CAMPOS2D_HXX
#define T8FLOW_CAMPOS2D_HXX

#include <cstddef>
#include <vector>

#include "malha3d.hxx"
#include "mls3d.hxx"

namespace t8flow3d {

//! Propriedade ESCALAR nos centros das celulas: um valor por folha.
//
// Indexada pelo id GLOBAL da folha, que e' o mesmo id que o solver usa -- de
// proposito: a tentacao de ter um indice proprio para os campos e' exatamente
// como os dois dominios do VOF passaram anos trocando indices entre si.
class CampoEscalar {
 public:
  explicit CampoEscalar(const Malha3 &m);

  //! Valor na folha de id global `id`.  Lanca se o id nao for local.
  double &em(long id);
  double em(long id) const;

  //! Valor na celula que CONTEM o ponto, pela busca geometrica da malha.
  //! Lanca `std::out_of_range` se o ponto estiver fora do dominio.
  double em_ponto(const Ponto3 &p) const;

  long tamanho() const { return static_cast<long>(v_.size()); }
  void preenche(double valor);

 private:
  const Malha3 *m_;
  long primeiro_;
  std::vector<double> v_;
};

//! Propriedade TENSORIAL nos centros: 2x2 componentes por folha.
//
// Uma componente por vez, como o `dpDu[i][j]` do HiGFlow -- e nao um `std::array`
// por celula, porque os lacos do metodo percorrem UMA componente de todas as
// celulas, nao todas as componentes de uma celula.  Layout component-major, que e'
// o que esse padrao de acesso pede.
class CampoTensor {
 public:
  explicit CampoTensor(const Malha3 &m);

  //! Componente (`lin`, `col`) na folha de id global `id`.
  double &em(long id, int lin, int col);
  double em(long id, int lin, int col) const;

  //! A componente inteira, como campo escalar nos centros -- para o laco do
  //! metodo poder percorrer uma componente sem saber do resto.
  CampoEscalar &componente(int lin, int col);
  const CampoEscalar &componente(int lin, int col) const;

  void preenche(double valor);

 private:
  // NOVE componentes em tres dimensoes, em ordem de linha: (0,0) (0,1) (0,2)
  // (1,0) ... (2,2).  Em duas eram quatro.
  std::vector<CampoEscalar> comp_;
};

//! A numeracao das FACETAS, uma lista por direcao.
//
// Produzida uma vez, na construcao, e consultada em O(1).  A ordem dos ids segue
// o percurso das folhas (a curva de preenchimento): para cada folha, primeiro a
// faceta "menos" daquela direcao e depois a "mais", quando esta for contorno.
class DominioFacetas {
 public:
  explicit DominioFacetas(const Malha3 &m);

  //! Facetas da direcao `d` (0 = x, 1 = y, 2 = z).
  long num_facetas(int d) const;

  //! Id da faceta da folha `(b, i)` no lado `lado` (0 = menos, 1 = mais) da
  //! direcao `d`, ou **-1 se esta folha nao for dona dela**.
  //
  // Devolver -1 e' informacao, nao falha: e' assim que o chamador sabe que a
  // faceta "mais" dele pertence ao vizinho, e o laco que soma fluxos tem de
  // consultar as DUAS pontas para nao contar faceta duas vezes.
  long id_faceta(long b, long i, int d, int lado) const;

  //! Centro geometrico da faceta `id` da direcao `d`.
  Ponto3 centro(int d, long id) const;

  //! A MEDIDA da faceta: em tres dimensoes uma AREA, o produto dos dois lados
  //! transversais da celula dona.
  //
  // Em duas dimensoes era um COMPRIMENTO, um lado so'.  Traduzir aquilo
  // literalmente daria um fluxo com a dimensao errada -- e um fluxo com
  // dimensao errada ainda converge, para a equacao errada.
  double medida(int d, long id) const;

  const Malha3 &malha() const { return *m_; }

 private:
  const Malha3 *m_;
  // por direcao, indexado pelo indice LOCAL da folha; -1 quando nao e' dona
  std::vector<long> menos_[3];
  std::vector<long> mais_[3];
  std::vector<Ponto3> centros_[3];
  std::vector<double> medidas_[3];
};

//! Propriedade VETORIAL nas facetas: uma componente por direcao, cada uma com um
//! valor por faceta daquela direcao.
//
// A componente `d` vive nas facetas da direcao `d` -- a componente x da
// velocidade nas facetas de normal x, e so' nelas.  Nao ha' "vetor completo numa
// faceta": isso e' o arranjo deslocado, e misturar as duas coisas e' como o
// arranjo colocado volta por descuido.
class CampoVetorial {
 public:
  explicit CampoVetorial(const DominioFacetas &fd);

  //! Componente `d` na faceta `id` daquela direcao.
  double &em(int d, long id);
  double em(int d, long id) const;

  void preenche(double valor);
  long tamanho(int d) const { return static_cast<long>(v_[d].size()); }

 private:
  const DominioFacetas *fd_;
  std::vector<double> v_[3];
};

}  // namespace t8flow3d

#endif
