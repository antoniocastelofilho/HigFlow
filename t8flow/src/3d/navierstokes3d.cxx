#include "navierstokes3d.hxx"

#include <cmath>
#include <stdexcept>
#include <string>

namespace t8flow3d {

NavierStokes::NavierStokes(const Malha3 &m, const DominioFacetas &fd,
                           const ContornoPorRegiao &bc_u0,
                           const ContornoPorRegiao &bc_u1,
                           const ContornoPorRegiao &bc_u2,
                           const ContornoPorRegiao &bc_p)
    : m_(&m), fd_(&fd), u_(fd), p_(m)
{
  bcu_[0] = &bc_u0;
  bcu_[1] = &bc_u1;
  bcu_[2] = &bc_u2;
  bcp_ = &bc_p;

  // AS DUAS RECUSAS, e nenhuma produz numero.
  //
  // Recusar na CONSTRUCAO, e nao no primeiro passo: quem monta uma malha
  // adaptada e so' descobre no passo mil que ela nao serve ja' gastou a
  // montagem, e pior, pode ter visto numeros no caminho.
  if (m.num_processos() > 1)
    throw std::invalid_argument(
        "NavierStokes3D: este solver e' SERIAL (" +
        std::to_string(m.num_processos()) +
        " processos).  Faltam a troca de franja dos campos -- o estencil le' "
        "celula e faceta do vizinho, e os campos cobrem so' os ids locais -- e "
        "o percurso por blocos, que supoe todas as arvores locais.  Com a "
        "franja ausente o estencil leria lixo perto da fronteira de particao e "
        "a corrida daria numero");

  const long trans = m.num_folhas_em_transicao();
  if (trans > 0)
    throw std::invalid_argument(
        "NavierStokes3D: a malha tem " + std::to_string(trans) +
        " folhas em transicao 2:1, e este solver ainda nao a trata.  Em tres "
        "dimensoes a face da celula grossa e' coberta por QUATRO facetas "
        "finas, nao duas, e a quina onde a viscosidade se interpola toca OITO "
        "celulas, nao quatro -- as regras de 2D aplicadas aqui produzem "
        "divergencia que nao fecha, e o sintoma aparece como o solver de "
        "pressao nao convergindo, longe da causa");

  // O MAPA de id global para `(bloco, folha)`.  O balanco de fluxos precisa
  // dele porque a faceta de um lado pode pertencer ao VIZINHO, e a' celula so'
  // sobra pedi-la a ele -- ver `divergencia_em`.
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      bi_por_id_[m.id_global(b, i)] = {b, i};

  u_.preenche(0.0);
  p_.preenche(0.0);
}

long NavierStokes::faceta_da_face(long b, long i, int d, int lado) const
{
  // A FACETA DE UMA FACE, VENHA ELA DE QUEM VIER.
  //
  // `DominioFacetas::id_faceta` devolve -1 quando ESTA folha nao e' dona da
  // faceta -- o que e' informacao e nao falha: e' assim que o chamador sabe
  // que ela pertence ao vizinho.  Mas o balanco de fluxos precisa das SEIS
  // faces da celula, donas ou nao: descartar as que nao sao dela deixa metade
  // do balanco de fora.
  //
  // MEDIDO AO ESQUECER DISTO: campo de divergencia 3 devolvia 18, e campo
  // solenoidal devolvia 10 em vez de zero.  Nao e' erro pequeno -- e' o
  // operador somando um subconjunto arbitrario de faces --, mas num caso de
  // ordem seria facil confundir com malha grossa.
  const long f = fd_->id_faceta(b, i, d, lado);
  if (f >= 0) return f;

  // Nao e' dona: a faceta e' a do OUTRO lado da celula vizinha por esta face.
  // A convencao de face e' a do t8code -- 0/1 em x, 2/3 em y, 4/5 em z.
  const std::vector<long> vz = m_->ids_vizinhos(b, i, 2 * d + lado);
  if (vz.size() != 1) return -1;      // contorno, ou transicao (recusada)
  const auto it = bi_por_id_.find(vz[0]);
  if (it == bi_por_id_.end()) return -1;
  return fd_->id_faceta(it->second.first, it->second.second, d, 1 - lado);
}

double NavierStokes::divergencia_em(long b, long i) const
{
  // O FLUXO LIQUIDO PELAS SEIS FACES, DIVIDIDO PELO VOLUME.
  //
  // A MEDIDA DA FACETA E' UMA AREA em tres dimensoes, e o volume e' o produto
  // dos tres lados.  Em malha uniforme `area/volume` e' `1/h` tanto quanto
  // seria com medida de comprimento e "volume" de area -- isto e', trocar as
  // duas por engano daria o MESMO numero aqui, e so' uma malha de celulas
  // nao cubicas separa as duas contas.  Por isso o caso que verifica isto usa
  // blocos de lados diferentes.
  const Ponto3 hh = m_->lados_da_folha(b, i);
  const double vol = hh[0] * hh[1] * hh[2];

  double fluxo = 0.0;
  for (int d = 0; d < 3; d++)
    for (int lado = 0; lado <= 1; lado++) {
      // O SINAL E' O DA NORMAL EXTERIOR: a face "menos" aponta para fora no
      // sentido negativo de `d`.
      const double sinal = (lado == 0) ? -1.0 : 1.0;
      const long f = faceta_da_face(b, i, d, lado);
      if (f < 0) continue;            // face de contorno sem faceta propria
      fluxo += sinal * u_.em(d, f) * fd_->medida(d, f);
    }
  return fluxo / vol;
}

double NavierStokes::divergencia_maxima() const
{
  double maxv = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++)
      maxv = std::max(maxv, std::fabs(divergencia_em(b, i)));
  return maxv;
}

}  // namespace t8flow3d
