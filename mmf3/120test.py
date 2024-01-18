import numpy as np
import matplotlib.pyplot as plt 
import sys

'''
 12. [5+5]
 a) Odredite potencijalnu energiju za zadanu konzervativnu silu F(x) = −2x, uzimaju´ ci nultu razinu u
 ishodištu. Na drugoj strani ovog papira u danoj mreži, napravite geometrijsku skicu procjene potencijalne
 energije pomo´ cu metode prediktor-korektor za 1 korak od ∆x = 1. Numerirajte redom me¯ dukorake, a
 derivacije u potrebnim toˇ ckama prikažite crtanjem kratkog segmenta pravca odgovaraju´ ceg nagiba.
 Komentirajte kakva bi bila greška procjene da smo napravili 100 puta više 100 puta kra´ cih koraka do
 toˇ cke x = 1.
 b) Nekajezadana parcijalna diferencijalna jednadžba (pojednostavljeni oblik kabelske jednadžbe koja opisuje
 prijenos akcijskog potencijala)
 λ2∂2V (x,t)
 ∂x2 −τ ∂V(x,t)
 =V(x,t)−A
 ∂t
 pri ˇ cemu su λ,τ i A konstante. Izvedite i napišite na papir jednadžbu u obliku prilago¯ denom za numeriˇ cko
 rješavanje eksplicitnom metodom.
'''


