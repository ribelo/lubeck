(ns ribelo.lubeck.redp
  (:require
   [ribelo.visby.math :as math]
   [ribelo.lubeck.quant :as quant]
   [ribelo.haag :as h]
   [net.cgrand.xforms :as x]))

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
  (^doubles [^double frisk ^double risk ^double freq ret]
   (let [arr (h/seq->double-array ret)
         redp (quant/rolling-economic-drawndown freq arr)
         std  (quant/annualized-risk freq)
         sr   (quant/annualized-sharpe-ratio frisk freq arr)]
     (math/min 1.0
               (math/max 0.0 (* (/ (+ (/ sr std) 0.5)
                                   (- 1.0 (math/pow risk 2.0)))
                                (math/max 0.0 (/ (- risk redp)
                                                 (- 1.0 redp)))))))))

(comment
  (do (quick-bench (into [] (redp-allocation 0.0 0.3 254) data))
      (quick-bench (redp-allocation 0.0 0.3 254 data))))
