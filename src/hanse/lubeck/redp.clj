(ns hanse.lubeck.redp
  (:require
   [uncomplicate.fluokitten.core :as fk]
   [hanse.halle :as h]
   [hanse.lubeck :as lubeck]
   [hanse.rostock.math :as math]
   [hanse.rostock.stats :as stats]))

(set! *unchecked-math* :warn-on-boxed)

(comment
  (do (def data (vec (repeatedly 100000 #(/ (- 0.5 ^double (rand)) 10.0))))
      (def arr  (double-array data))
      (require '[criterium.core :refer [quick-bench]])
      (require '[taoensso.encore :as e])))

(defn single-allocation
  ^double [^double frisk ^double risk ^long freq close]
  (let [close' (h/seq->double-array close)
        ret    (lubeck/tick->ret close')
        redp   (lubeck/rolling-economic-drawndown freq close')
        std    (lubeck/annualized-risk freq ret)
        sr     (lubeck/annualized-sharpe-ratio frisk freq ret)]
    (math/min 1.0
              (math/max 0.0 (* (/ (+ (/ sr std) 0.5)
                                  (- 1.0 (math/pow risk 2.0)))
                               (math/max 0.0 (/ (- risk redp)
                                                (- 1.0 redp))))))))

(comment
  (single-allocation 0.0 0.3 5 data))

(defn- redp-stats [^double frisk ^double risk ^long freq ^doubles close]
  (let [close'  (h/seq->double-array close)
        ret     (lubeck/tick->ret close')
        std     (lubeck/annualized-risk freq ret)
        ann-ret (lubeck/ann-return-geometric freq ret)
        drift   (math/max 0.0 (+ (- ann-ret frisk) (/ (math/sq std) 2.0)))
        redp    (lubeck/rolling-economic-drawndown freq close')
        Y       (* (/ 1.0 (- 1.0 (math/sq risk)))
                   (/ (- risk redp)
                      (- 1.0 redp)))]
    [std drift Y]))

(comment
  (redp-stats 0.0 0.3 5 data))

(defn multiple-allocation
  ([frisk risk freq assets]
   (let [[volatility
          drift
          Y]            (apply map vector (fk/fmap (partial redp-stats frisk risk freq) assets))
         inv-volatility (fk/fmap #(math/pow % -1.0) volatility)
         matrix         (fk/fmap (comp math/round *) inv-volatility drift inv-volatility Y)
         sum            (fk/fold + matrix)]
     (if (<= ^double sum 1.0) matrix (fk/fmap #(/ ^double % ^double sum) matrix)))))
