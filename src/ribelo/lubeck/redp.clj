(ns ribelo.lubeck.redp
  (:require
   [primitive-math :as p]
   [ribelo.haag :as h]
   [ribelo.lubeck.quant :as quant]
   [ribelo.visby.math :as math]
   [uncomplicate.fluokitten.core :as fk]))

(set! *warn-on-reflection* true)

(comment
  (do (def data (vec (repeatedly 100000 #(/ (- 0.5 (rand)) 10.0))))
      (def arr  (double-array data))))

(defn single-allocation
  ^double [^double frisk ^double risk ^long freq close]
  (let [close' (h/seq->double-array close)
        ret    (quant/tick->ret close')
        redp   (quant/rolling-economic-drawndown freq close')
        std    (quant/annualized-risk freq ret)
        sr     (quant/annualized-sharpe-ratio frisk freq ret)]
    (math/min 1.0
              (math/max 0.0 (p/* (p/div (p/+ (p/div sr std) 0.5)
                                        (p/- 1.0 (math/pow risk 2.0)))
                                 (math/max 0.0 (p/div (p/- risk redp)
                                                      (p/- 1.0 redp))))))))

(comment
  (quick-bench (single-allocation 0.0 0.3 5 data)))

(defn- redp-stats [^double frisk ^double risk ^long freq close]
  (let [ret     (quant/tick->ret close)
        std     (quant/annualized-risk freq ret)
        ann-ret (quant/ann-return-geometric freq ret)
        drift   (math/max 0.0 (p/+ (p/- ann-ret frisk) (p/div (math/sq std) 2.0)))
        reddp   (quant/rolling-economic-drawndown freq close)
        Y       (p/*
                 (p/div 1.0 (p/- 1.0 (math/sq risk)))
                 (p/div (p/- risk reddp)
                        (p/- 1.0 reddp)))]
    [std drift Y]))

(defn multiple-allocation [frisk risk freq assets]
  (let [symbols        (into [] (map #(:symbol (first %))) assets)
        [volatility
         drift
         Y]            (apply map vector (fk/fmap (partial redp-stats frisk risk freq) assets))
        inv-volatility (fk/fmap #(math/pow % -1.0) volatility)
        matrix         (fk/fmap (comp math/round *) inv-volatility drift inv-volatility Y)
        sum            (double (fk/fold + matrix))
        proportions    (if (p/<= sum 1.0) matrix (fk/fmap #(/ % sum) matrix))]
    (fk/fmap
     (fn [symbol p]
       {:symbol     symbol
        :allocation p})
     symbols proportions)))
