#= Date Wikipedia (abundente) + IAEA (sectiuni eficace la 0.25 eV) 
https://www-nds.iaea.org/ngatlas2/ =#

abundente_izotopice = [1.25, 0.89, 12.47, 12.80, 24.11, 12.23, 28.75, 7.51]
abundente_izotopice = abundente_izotopice./100 # Transformam in procente

sectiuni_eficace = [0.30, 0.34, 3.50, 6.80, 0.60, 18000.0, 0.05, 0.015] #in barni

σCd = sum(abundente_izotopice .* sectiuni_eficace)