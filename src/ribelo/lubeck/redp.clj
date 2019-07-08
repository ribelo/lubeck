(ns ribelo.lubeck.redp
  (:require
   [ribelo.visby.math :as math]
   [ribelo.lubeck.quant :as quant]
   [ribelo.haag :as h]
   [net.cgrand.xforms :as x]))

(comment
  (do (def data (vec (repeatedly 100000 #(/ (- 0.5 (rand)) 10.0))))
      (def arr  (double-array data))))

(defn single-allocation
  ([frisk risk freq]
   (comp
    (x/transjuxt [(quant/rolling-economic-drawndown freq)
                  (quant/annualized-risk freq)
                  (comp (quant/tick->ret) (quant/annualized-sharpe-ratio frisk))])
    (x/reduce
     (fn
       ([] [0.0 0.0 0.0])
       ([[redp std sharpe-ratio]]
        (math/min 1.0
                  (math/max 0.0 (* (/ (+ (/ sharpe-ratio std) 0.5)
                                      (- 1.0 (math/pow risk 2.0)))
                                   (math/max 0.0 (/ (- risk redp)
                                                    (- 1.0 redp)))))))
       ([acc coll] coll)))))
  (^double [^double frisk ^double risk ^long freq close]
   (let [close'  (h/seq->double-array close)
         ret  (quant/tick->ret close')
         redp (quant/rolling-economic-drawndown freq close')
         std  (quant/annualized-risk freq ret)
         sr   (quant/annualized-sharpe-ratio frisk freq ret)]
     (math/min 1.0
               (math/max 0.0 (* (/ (+ (/ sr std) 0.5)
                                   (- 1.0 (math/pow risk 2.0)))
                                (math/max 0.0 (/ (- risk redp)
                                                 (- 1.0 redp)))))))))

(comment
  (do (quick-bench (into [] (single-allocation 0.0 0.3 254) data))
      (quick-bench (single-allocation 0.0 0.3 100000 data))))
