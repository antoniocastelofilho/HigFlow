Modelo Integral KBKZ PSM ou UCM  Integral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OBS:Para salvar pressao, tensores e velocidades para visualizacao fora do paraview 
crie as pastas Dados/Pressao  Dados/Tensores e Dados/Velocidade
e habilite as rotinas de impressao no arquivo ns-exemple-2d.c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     arquivo de entrada "nome.viscintcontr"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
São 4 parâmetros necessários:

contr.model       =======> [O] Modelo Integral
                  =======> [1] Integral Fracionário

contr.model_H     =======> [O] Modelo PSM  
                  =======> [1] Modelo UCM (integral)

contr.discrtype   =======> [O] Euler Explícito
 

convecdiscrtype   =======> [O] Esquema Central
                  =======> [1] Esquema Cubista

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     arquivo de entrada "nome.viscintpar"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
São 25 parâmetros necessários:

ns->ed.im.par.De             =======> Deborah;

ns->ed.im.par.alpha          =======> Função H;
ns->edim.par.beta            =======> Função H;

ns->ed.im.par.scorte         =======> intervalo de tempo (passado) 
                                      deve ser da ordem do tempo de relaxamento 
                                      lambda_ref(NUNCA USE ZERO!!!!!)

ns->ed.im.par.rho            =======> Densidade

ns->ed.im.par.v_ref          =======> Escala Velocidade
ns->ed.im.par.l_ref          =======> Lambda_ref

% parâmetros específicos de cada fluido

ns->ed.im.par.M              =======> Quantidade de módulos

Para k=0:7  
ns->ed.im.par.a[k]           =======> Modulo relaxamento
ns->ed.im.par.lambda[k]      =======> Tempo relaxamento
  
OBS: se M<7 então para k>=M, a[k]=0 e lambda[k]=1;
     Os valores de a[k] e lambda[s] são DIMENSIONAIS 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
