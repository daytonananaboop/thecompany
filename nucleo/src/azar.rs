pub struct Azar {
    estado: u64,
}

impl Azar {
    pub fn nuevo(semilla: u64) -> Azar {
        Azar { estado: semilla }
    }

    pub fn siguiente(&mut self) -> u64 {
        self.estado = self.estado.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.estado;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn decimal(&mut self) -> f64 {
        (self.siguiente() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub fn entero(&mut self, hasta: usize) -> usize {
        if hasta == 0 {
            return 0;
        }
        (self.decimal() * hasta as f64) as usize % hasta
    }

    pub fn elegir_con_pesos(&mut self, pesos: &[f64]) -> usize {
        let total: f64 = pesos.iter().sum();
        if pesos.is_empty() || total <= 0.0 {
            return self.entero(pesos.len());
        }
        let mut objetivo = self.decimal() * total;
        for (indice, peso) in pesos.iter().enumerate() {
            if objetivo < *peso {
                return indice;
            }
            objetivo -= peso;
        }
        pesos.len() - 1
    }

    pub fn normal(&mut self) -> f64 {
        let u = self.decimal().max(f64::MIN_POSITIVE);
        let v = self.decimal();
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
    }
}
