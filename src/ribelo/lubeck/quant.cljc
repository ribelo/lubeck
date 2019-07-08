(ns ribelo.lubeck.quant
  (:require
   [clojure.spec.alpha :as s]
   [net.cgrand.xforms :as x]
   [net.cgrand.xforms.rfs :as rf]
   #?(:clj [uncomplicate.fluokitten.core :as fk])
   #?(:clj [uncomplicate.fluokitten.jvm])
   ;; [criterium.core :refer [quick-bench]]
   [ribelo.haag :as h]
   [ribelo.visby.math :as math]
   [ribelo.visby.emath :as emath]
   [ribelo.visby.stats :as stats])
  (:import (org.apache.commons.math3.stat StatUtils)))

(comment
  (do (def data (vec (repeatedly 100000 #(/ (- 0.5 (rand)) 10.0))))
      (def arr  (double-array data))))

(defn ann-return-geometric
  ([^double freq]
   (comp
    (x/reduce
     (fn
       ([] [])
       ([[x n]] (- (math/pow x (/ freq n)) 1.0))
       ([acc [x n]] [x n])))))
  ([^double freq ret]
   (let [arr    (h/seq->double-array ret)
         n      (alength arr)
         return (StatUtils/product (emath/add 1.0 arr))]
     (- (math/pow return (/ freq n)) 1.0))))

(comment
  (do (quick-bench (into [] (ann-return-geometric 254) data))
      (quick-bench (ann-return-geometric 254.0 data))))

(defn ann-return-simple
  ([^double freq]
   (comp
    (stats/mean)
    (emath/mul freq)))
  ([^double freq ret]
   (let [arr (h/seq->double-array ret)
         mean (stats/mean arr)]
     (* mean freq))))

(comment
  (do (quick-bench (into [] (ann-return-simple 254.0) data))
      (quick-bench (ann-return-simple 254.0 data))))

(defn annualized-return
  "Average annualized returns over a period, convenient when comparing returns.
  It can be an Arithmetic or Geometric (default) average return: if compounded with itself the
  geometric average will be equal to the cumulative return"
  ([^double freq mode]
   (case mode
     :geometric (ann-return-geometric freq)
     :simple    (ann-return-simple freq)))
  ([^double freq mode ret]
   (case mode
     :geometric (ann-return-geometric freq ret)
     :simple    (ann-return-simple freq ret))))

(defn active-return ;; TODO
  "Asset/Portfolio annualized return minus Benchmark annualized return"
  ([xs freq mode]
   (let [xs' (sequence (annualized-return freq mode) xs)]
     (comp
      (annualized-return freq mode)
      (emath/sub xs')))))

(defn annualized-risk
  "Annualized standard deviation of asset/portfolio returns"
  ([^double freq]
   (comp
    (stats/std)
    (emath/mul (math/sqrt freq))))
  ([^double freq ret]
   (let [arr (h/seq->double-array ret)]
     (* (stats/std arr) (math/sqrt freq)))))

(comment
  (do (quick-bench (into [] (annualized-risk 254.0) data))
      (quick-bench (annualized-risk 254.0 data))))

(defn sharpe-ratio
  "Sharpe Ratio.Compute Sharpe ratio for an collection XS of values (daily, weekly, etc) and
   a free-risk rate. Annual free-risk must be divided to match the right timeframe."
  ([^double frisk]
   (comp
    (x/transjuxt [(stats/mean) (stats/std-s)])
    (x/reduce
     (fn
       ([] [0.0 0.0])
       ([[mean std]] (/ (- mean frisk) std))
       ([acc coll] coll)))))
  ([^double frisk ^double freq ret]
   (let [arr (h/seq->double-array ret)
         mean (stats/mean arr)
         std (stats/std-s arr)]
     (/ (- mean frisk) std))))

(comment
  (do (quick-bench (into [] (sharpe-ratio 0.0 254.0) data))
      (quick-bench (sharpe-ratio 0.0 254.0 data))))

(defn annualized-sharpe-ratio
  ([^double frisk]
   (comp
    (x/transjuxt [(annualized-return) (annualized-risk)])
    (x/reduce
     (fn
       ([] [0.0 0.0])
       ([[mean std]] (/ (- mean frisk) std))
       ([acc coll] coll)))))
  (^double [frisk ^doubles ret]
   (let [ann-ret (annualized-return ret)
         std     (annualized-risk ret)]
     (/ (- ann-ret frisk) std))))

(defn adjusted-sharpe-ratio ;;TODO
  "Sharpe Ratio adjusted for skewness and kurtosis with a penalty factor
   for negative skewness and excess kurtosis."
  ([^double frisk]
   (comp
    (x/transjuxt [(sharpe-ratio frisk) stats/skewness stats/kurtosis])
    (x/reduce
     (fn
       ([] (transient []))
       ([[sr sk ku]] (* sr (- (+ 1 (* (/ sk 6) sr))
                              (* (/ (- ku 3) 24) (math/sqrt sr)))))
       ([acc [sr sk ku]] (-> acc (conj! sr) (conj! sk) (conj! ku)))))))
  ([] (adjusted-sharpe-ratio 0.0)))

(defn annualized-adjusted-sharpe-ratio ;;TODO
  "Sharpe Ratio adjusted for skewness and kurtosis with a penalty factor
   for negative skewness and excess kurtosis."
  ([^double frisk ^long freq mode]
   (comp
    (x/transjuxt [stats/skewness stats/kurtosis
                  (annualized-return freq mode) (annualized-risk freq)])
    (x/reduce
     (fn
       ([] (transient []))
       ([[sk ku annret annrisk]]
        (let [sr (/ (- (/ (Math/round ^double (* 10000 annret)) 10000) frisk)
                    (/ (Math/round ^double (* 10000 annrisk)) 10000))]
          (* sr (- (+ 1 (* (/ sk 6) sr))
                   (* (- ku 3) (math/sqrt sr))))))
       ([acc [sk ku annret annrisk]]
        (-> acc (conj! sk) (conj! ku) (conj! annret) (conj! annrisk)))))))
  ([^double frisk ^long freq]
   (annualized-adjusted-sharpe-ratio frisk freq :geometric))
  ([^double frisk]
   (annualized-adjusted-sharpe-ratio frisk 252 :geometric))
  ([]
   (annualized-adjusted-sharpe-ratio 0.0 252 :geometric)))

(defn downside-risk
  "Downside Risk or Semi-Standard Deviation.
   Measures the variability of underperformance below a minimum target rate"
  ([^double mar]
   (let [fnx (comp
              (emath/sub mar)
              (emath/min 0.0)
              (emath/pow 2.0)
              (x/into []))
         fny (fn
               ([] [])
               ([[xs n]]
                (let [sum (first (into [] (comp (emath/div (double n)) (x/reduce +)) xs))]
                  (math/sqrt ^double sum)))
               ([acc [xs count]] (-> acc (conj xs) (conj count))))]
     (comp
      (x/transjuxt [fnx x/count])
      (x/reduce fny))))
  ([^double mar ret]
   (let [arr (h/seq->double-array ret)
         n   (alength ^doubles arr)]
     (loop [i 0 sum 0.0]
       (if (< i n)
         (recur (inc i) (+ sum
                           (/ (math/sq (math/min 0.0 (- (aget ^doubles arr i) mar)))
                              n)))
         (math/sqrt sum))))))

(comment
  (do (quick-bench (into [] (downside-risk 0.0) data))
      (quick-bench (downside-risk 0.0 data))))

(defn sortino-ratio
  "Sortino ratio"
  ([^double frisk ^double mar]
   (comp
    (x/transjuxt [(downside-risk mar) (stats/mean)])
    (x/reduce
     (fn
       ([] (transient []))
       ([[mean dr]] (/ (- mean frisk) dr))
       ([acc [dr mean]] (-> acc (conj! mean) (conj! dr)))))))
  ([^double frisk ^double mar ret]
   (let [arr (h/seq->double-array ret)
         dr (downside-risk mar arr)
         mean (stats/mean arr)]
     (/ (- mean frisk) dr))))

(comment
  (do (quick-bench (into [] (sortino-ratio 0.0 0.0) data))
      (quick-bench (sortino-ratio 0.0 0.0 data))))

(defn drawdown
  "Drawdowon from Peak. Any continuous losing return period."
  ([]
   (comp
    (emath/add 1.0)
    (emath/cumprod)
    (fn [xf]
      (let [val (volatile! 1.0)]
        (fn
          ([] (xf))
          ([acc] (xf acc))
          ([acc x]
           (let [val' (math/max x @val)]
             (vreset! val val')
             (xf acc (/ (- val' x) val')))))))))
  ([close-or-ret]
   (let [arr (h/seq->double-array close-or-ret)
         n   (alength ^doubles arr)
         r   (double-array n)]
     (loop [i 0 s 1.0 mx 1.0]
       (if (< i n)
         (let [v (* s (+ 1.0 (aget ^doubles arr i)))
               mx' (math/max v mx)
               dr (/ (- mx' v) mx')]
           (aset r i dr)
           (recur (inc i) v mx'))
         r)))))

(comment
  (do (quick-bench (into [] (drawdown) data))
      (quick-bench (drawdown data))))

(defn continuous-drawdown
  ([]
   (comp
    (partition-by neg?)
    (filter #(every? neg? %))
    (fn [rf]
      (fn
        ([] (rf))
        ([acc] (rf acc))
        ([acc x] (rf acc (- 1 (transduce (map #(+ 1 %)) * x))))))))
  ([close-or-ret]
   (let [arr (h/seq->double-array close-or-ret)
         n   (alength ^doubles arr)
         dq   (java.util.ArrayDeque.)]
     (loop [i 0 s 1.0]
       (if (< i n)
         (let [v (aget ^doubles arr i)]
           (cond
            (and (zero? i) (< v 0.0))
            (recur (inc i) (+ 1.0 v) )
            (< 0.0 i)
            (if (< v 0.0)
              (recur (inc i) (* s (+ 1.0 v)))
              (let [dd (- 1.0 s)]
                (when-not (zero? dd)
                  (.add dq dd))
                (recur (inc i) 1.0)))))
         (let [r (double-array (.toArray dq))]
           (.clear dq)
           r))))))

(comment
  (do (quick-bench (into [] (continuous-drawdown) data))
      (quick-bench (continuous-drawdown data))))

(defn average-drawdown
  ([]
   (comp (continuous-drawdown)
         (stats/mean)))
  ([close-or-ret]
   (let [arr (h/seq->double-array close-or-ret)]
     (->> (continuous-drawdown arr)
          (stats/mean)))))

(comment
  (do (quick-bench (into [] (average-drawdown) data))
      (quick-bench (average-drawdown data))))

(defn maximum-drawdown
  ([] (comp (continuous-drawdown)
            (stats/max)))
  ([close-or-ret]
   (let [arr (h/seq->double-array close-or-ret)]
     (->> (continuous-drawdown arr)
          (stats/max)))))

(comment
  (do (quick-bench (into [] (continuous-drawdown) data))
      (quick-bench (continuous-drawdown data))))

(def rate-of-return ;;TODO
  "Simple rate of return calculated from the last and the first value of
  an array of numbers."
  (comp (emath/add 1.0)
        (x/reductions *)
        (x/transjuxt [(x/reduce rf/some)
                      x/last])
        (x/reduce
         (fn
           ([] (transient []))
           ([[f l]] (- (/ l f) 1.0))
           ([acc [f l]] (-> acc (conj! f) (conj! l)))))))

(defn cagr ;;TODO
  "Compound annual growth rate"
  ([n]
   (comp rate-of-return
         (emath/add 1.0)
         (emath/pow (/ 1.0 n))
         (emath/sub 1.0)))
  ([] rate-of-return))

(defn calmar-ratio ;;TODO
  "A risk-adjusted measure like Sharpe ratio that uses maximum drawdown instead of
  standard deviation for risk."
  ([^double frisk ^long freq]
   (comp (x/transjuxt [(annualized-return freq :geometric) maximum-drawdown])
         (x/reduce
          (fn
            ([] (transient []))
            ([[annret maxdd]] (/ (- annret frisk) maxdd))
            ([acc [annret maxdd]] (-> acc (conj! annret) (conj! maxdd)))))))
  ([frisk]
   (calmar-ratio frisk 252))
  ([]
   (calmar-ratio 0.0 252)))

(defn downside-potential ;;TODO
  "Downside potential"
  ([mar]
   (comp (x/transjuxt [(comp (emath/mul -1.0) (emath/add mar) (emath/max 0.0) (x/into [])) x/count])
         (x/reduce
          (fn
            ([] (transient []))
            ([[xs count]] (transduce (emath/div count) + xs))
            ([acc [xs count]] (-> acc (conj! xs) (conj! count)))))))
  ([] (downside-potential 0.0)))

(defn burke-ratio ;;TODO
  "A risk-adjusted measure with free risk and drawdowns.
   For the 'simple' mode the excess return over free risk is divided by the square root of
   the sum of the square of the drawdowns. For the 'modified' mode the Burke Ratio is multiplied
   by the square root of the number of datas."
  ([frisk freq mode]
   (comp
    (x/transjuxt [(annualized-return freq)
                  (comp continuous-drawdown
                        (emath/pow 2)
                        (x/reduce +)
                        emath/sqrt)
                  x/count])
    (x/reduce
     (fn
       ([] (transient []))
       ([[annret dd c]] (case mode
                          :simple (/ (- annret frisk) dd)
                          :modified (* (/ (- annret frisk) dd)
                                       (math/sqrt c))))
       ([acc [annret dd c]] (-> acc
                                (conj! annret)
                                (conj! dd)
                                (conj! c)))))))
  ([frisk freq]
   (burke-ratio frisk freq :simple))
  ([frisk]
   (burke-ratio frisk 252 :simple))
  ([]
   (burke-ratio 0.0 252 :simple)))

(def ulcer-index ;;TODO
  "Ulcer Index of Peter G. Martin (1987). The impact of long, deep drawdowns will have significant
  impact because the underperformance since the last peak is squared."
  (comp
   (x/transjuxt [(comp drawdown
                       (emath/pow 2)
                       (x/reduce +))
                 x/count])
   (x/reduce
    (fn
      ([] (transient []))
      ([[dd c]] (math/sqrt (/ dd c)))
      ([acc [dd c]] (-> acc (conj! dd) (conj! c)))))))

(defn martin-ratio ;;TODO
  "A risk-adjusted measure with free risk and Ulcer index.
   Martin Ratio = (Portfolio Return - RiskFree) / Ulcer Index
   Mode: :return, :geometric (default: :return)"
  ([frisk freq]
   (comp
    (x/transjuxt [(annualized-return freq)
                  ulcer-index])
    (x/reduce
     (fn
       ([] (transient []))
       ([[annret u]] (/ (- annret frisk) u))
       ([acc [annret u]] (-> acc (conj! annret) (conj! u)))))))
  ([frisk]
   (martin-ratio frisk 252))
  ([]
   (martin-ratio 0.0 252)))

(def hurst-index ;;TODO
  "It's a useful statistic for detecting if a time series is mean reverting (anti-persistent), totally random or persistent.
   A value in the range [0.5) indicates mean-reverting (anti-persistent)
   A value of 0.5 indicate a random walk
   A value H in the range (0.5,1] indicates momentum (persistent)"
  (comp
   (x/transjuxt [(comp emath/cumdev stats/max)
                 (comp emath/cumdev stats/min)
                 stats/std
                 x/count])
   (x/reduce
    (fn
      ([] (transient []))
      ([[mx mn std n]]
       (let [rs (/ (- mx mn) std)]
         (/ (math/log rs) (math/log n))))
      ([acc [mx mn std n]] (-> acc
                               (conj! mx) (conj! mn)
                               (conj! std) (conj! n)))))))

(defn info-ratio ;;TODO
  "Information Ratio"
  [xs]
  (comp
   (x/transjuxt [(comp (emath/sub xs) stats/std)
                 (comp (emath/sub xs) stats/mean)])
   (x/reduce
    (fn
      ([] (transient []))
      ([[std mean]] (/ mean std))
      ([acc [std mean]] [std mean])))))

(defn jensen-alpha ;;TODO
  "Ex-post alpha calculated with regression line.
  Free-risk is the avereage free-risk for the timeframe selected."
  [xs frisk]
  (comp
   (x/transjuxt [(comp (emath/sub xs) stats/std)
                 (comp (emath/sub xs) stats/mean)])
   (x/reduce
    (fn
      ([] (transient []))
      ([[std mean]]
       (/ mean std))
      ([acc [std mean]] (-> acc (conj! std) (conj! mean)))))))

(defn modigliani ;;TODO
  "Modigliani index for risk-adjusted return"
  [ys frisk]
  (let [[stdb] (sequence stats/std ys)]
    (comp
     (x/transjuxt [stats/mean
                   (sharpe-ratio frisk)
                   stats/std])
     (x/reduce
      (fn
        ([] [0.0 0.0 0.0])
        ([[mean sharpe std]]
         (+ mean (* sharpe (- stdb std))))
        ([acc coll] coll))))))

(defn rolling-economic-drawndown
  ([^long freq]
   (comp
    (x/take-last freq)
    (x/transjuxt [(x/take-last 1)
                  (stats/max)])
    (map (fn [[last max]] (double (- 1.0 (/ last max)))))))
  ([^long freq close]
   (let [arr (-> (h/seq->double-array close) (h/take-last freq))
         max (stats/max ^doubles arr)]
     (- 1.0 (/ ^double (h/last arr) max)))))

(comment
  (do (quick-bench (into [] (rolling-economic-drawndown 254) data))
      (quick-bench (rolling-economic-drawndown 254 data))))

(defn tick->ret
  "Convert a value series to a return series"
  ([]
   (comp (x/partition 2 1)
         (fn [rf]
           (fn
             ([] (rf))
             ([acc] (rf acc))
             ([acc [x y]] (rf acc (- (/ y x) 1.0)))))))
  ([close]
   (let [arr (h/seq->double-array close)
         n   (alength ^doubles arr)
         nn  (dec n)
         r   (double-array (dec n))]
     (loop [i 0]
       (when (< i nn)
         (let [x (aget ^doubles arr i)
               y (aget ^doubles arr (inc i))
               ret (- (/ y x) 1.0)]
           (aset ^doubles r i ret)
           (recur (inc i)))))
     r)))

(comment
  (do (quick-bench (into [] (tick->ret) data))
      (quick-bench (tick->ret data))))
