#lang racket
(require 2htdp/universe racket/draw racket/generator images/logos)
;----------- const --------------------
(define log2 (log 2))
(define (sqr/sc x) (/ (sqr x) (scale)))
;----------- independent parameters -----------
(define zoom (make-parameter 1/2))
(define limit (make-parameter 100))
(define bailout (make-parameter 100))
(define center (make-parameter 0))
;-------------------- others --------------------------
(define radius (make-derived-parameter bailout sqr sqrt))
(define view (make-parameter 600+450i))
(define scale (make-parameter 225))
(define (set-size! x)
  (scale (* 1/2 (min (real-part x) (imag-part x))))
  (view x))

(define sac-stripe (make-parameter 5))
(define sac-trunc (make-parameter 1))
;----------------------- fractal types ---------------------------
(define (mandelbrot-gen c)
  (let ([i (/ 2 (scale))])
    (generator ()
               (let loop ([z (exact->inexact c)] [x2 (exact->inexact (sqr/sc (real-part c)))]
                          [y2 (exact->inexact (sqr/sc (imag-part c)))])
               (yield z x2 y2)
               (let ([t (+ c (make-rectangular (- x2 y2) (* i (real-part z) (imag-part z))))])
                 (loop t (sqr/sc (real-part t)) (sqr/sc (imag-part t))))))))
(define (bship-gen c)
  (let ([i (/ 2 (scale))])
    (generator ()
               (let loop ([z (exact->inexact c)] [x2 (exact->inexact (sqr/sc (real-part c)))]
                          [y2 (exact->inexact (sqr/sc (imag-part c)))])
               (yield z x2 y2)
               (let ([t (+ c (make-rectangular (- x2 y2) (* -1 (abs (* i (real-part z) (imag-part z))))))])
                 (loop t (sqr/sc (real-part t)) (sqr/sc (imag-part t))))))))
(define (julia-gen p)
  (let ([j (exact->inexact (* p (/ (scale) (zoom))))] [i (/ 2 (scale))])
  (λ (c) (generator ()
                    (let loop ([z (exact->inexact c)] [x2 (exact->inexact (sqr/sc (real-part c)))]
                               [y2 (exact->inexact (sqr/sc (imag-part c)))])
                      (yield z x2 y2)
                      (let ([t (+ j (make-rectangular (- x2 y2) (* i (real-part z) (imag-part z))))])
                        (loop t (sqr/sc (real-part t)) (sqr/sc (imag-part t)))))))))
(define (bship-test c)
  (or (and (zero? (imag-part c))
           (<= (real-part c) (/ (scale) 4))
           (>= (real-part c) (* -2 (scale))))
      (and (positive? (imag-part c))
           (<= (real-part c) (/ (scale) 4))
           (let ([ys (sqr (imag-part c))])
             (or (< (/ (+ (sqr (+ (real-part c) (scale))) ys) (scale)) (/ (scale) 16))
                 (let* ([t (- (real-part c) (/ (scale) 4))]
                        [q (/ (+ (sqr t) ys) (scale))])
                   (< (* q (+ q t)) (/ ys 4))))))))
(define (mandelbrot-test c)
  (or (and (zero? (imag-part c))
       (<= (real-part c) (/ (scale) 4))
       (>= (real-part c) (* -2 (scale))))
      (let ([ys (sqr (imag-part c))])
        (or (< (/ (+ (sqr (+ (real-part c) (scale))) ys) (scale)) (/ (scale) 16))
            (let* ([t (- (real-part c) (/ (scale) 4))]
                   [q (/ (+ (sqr t) ys) (scale))])
              (< (* q (+ q t)) (/ ys 4)))))))
(define-syntax (mandelbrot-escape stx)
  (syntax-case stx (dem)
    [(_ dem) #'(dem mandelbrot-gen mandelbrot-test 1 1 1)]
    [(_ esc args ...) #'(esc mandelbrot-gen mandelbrot-test 1 args ...)]))
(define-syntax (julia-escape stx)
  (syntax-case stx (dem)
    [(_ p dem) #'(dem (julia-gen p) (λ (x) #f) 0 1 0)]
    [(_ p esc args ...) #'(esc (julia-gen p) (λ (x) #f) 0 args ...)]))
(define-syntax (bship-escape stx)
  (syntax-case stx (dem)
    [(_ dem) #'(dem bship-gen bship-test 1 1 1)]
    [(_ esc args ...) #'(esc bship-gen bship-test 1 args ...)]))

;--------------------- escapes ----------------------------------
(define (normal gen test start)
  (λ (c) (if (test c) 0 (- 1 (/ (point-loop (n start) (gen z x2 y2) #f c n) (limit))))))
(define (potential gen test start)
  (let ([ls (log (scale))] [b (* (scale) (bailout))])
    (λ (c) (if (test c) 0 ((λ (p) (if (<= (cdr p) b) 0
                                      (expt (* 4 (/ (- (log (cdr p)) ls)
                                                    (arithmetic-shift 1 (add1 (car p)))))
                                            (/ 1 (+ 7 (integer-length (exact-floor (zoom))))))))
                           (point-loop (n start) (gen z x2 y2) #f c (cons n (+ x2 y2))))))))
(define (smooth-par b b2 lb ls)
  (λ (mag) (cond [(<= mag b) 1]
                 [(>= mag b2) 0]
                 [else (add1 (- lb (/ (log (- (log mag) ls)) log2)))])))
(define (smooth gen test start)
  (let ([par (smooth-par (* (scale) (bailout)) (* (scale) (sqr (bailout))) (/ (log (log (bailout))) log2) (log (scale)))])
    (λ (c) (if (test c) 0
               ((λ (p) (- 1 (/ (+ (car p) (par (cdr p))) (add1 (limit)))))
                (point-loop (n start) (gen z x2 y2) #f c (cons n (+ x2 y2))))))))
(define (dem-out ls b)
  (λ (d mag)
    (if (or (zero? d) (<= mag b)) 0
        (expt (* 4 (- (log mag) ls)
                 (sqrt (/ mag (* (scale) (+ (sqr (real-part d)) (sqr (imag-part d)))))))
              (/ 1 (+ 3 (integer-length (exact-floor (zoom)))))))))
(define (dem gen test start [init-d 1] [add 1])
  (let ([i (/ 2 (scale))] [out (dem-out (log (scale)) (* (scale) (bailout)))])
    (λ (c) (let ([d init-d] [next-d init-d])
             (if (test c) 0
                 ((λ (mag) (out d mag))
                  (point-loop (n start) (gen z x2 y2) #f c
                              (set! d next-d) (set! next-d (+ add (* i d z))) (+ x2 y2))))))))
(define (sac-t z)
  (* 1/2 (add1 (sin (* (sac-stripe) (angle z))))))
(define ((sac-out par) a-1 n z mag)
  (if (or (>= n (limit)) (< n (sac-trunc))) 0
      (let ([d (par mag)])
        (+ (* d (/ (+ a-1 (sac-t z)) (- n -1 (sac-trunc))))
           (if (> n (sac-trunc)) (* (- 1 d) (/ a-1 (- n (sac-trunc)))) 0)))))
(define (sac gen test start)
  (let* ([b (* (scale) (bailout))]
         [out (sac-out (smooth-par b (* (scale) (sqr (bailout)))
                                   (/ (log (log (bailout))) log2) (log (scale))))])
    (λ (c) (let ([a-1 0])
             (if (test c) 0
                 ((λ (lst) (apply out a-1 lst))
                  (point-loop (n start) (gen z x2 y2) #f c
                              (let ([mag (+ x2 y2)])
                                (unless (or (< n (sac-trunc)) (> mag b))
                                  (set! a-1 (+ a-1 (sac-t z))))
                                (list n z mag)))))))))
(define ((point-trap up) gen test start)
  (let ([p (* (/ (scale) (zoom)) up)])
    (λ (c) (let ([trap +inf.0] [x +inf.0] [t +inf.0])
             (point-loop (n start) (gen z x2 y2) #f c
                         (set! t (- z p))
                         (set! x (+ (sqr/sc (real-part t)) (sqr/sc (imag-part t))))
                         (when (< x trap) (set! trap x)))
             (sqrt (/ trap (scale)))))))
(define ((line-trap up) gen test start)
  (let ([p (* (/ (scale) (zoom)) up)])
    (λ (c) (let ([trap +inf.0] [x +inf.0])
             (point-loop (n start) (gen z x2 y2) #f c
                         (set! x (min (abs (- (real-part z) (real-part p))) (abs (- (imag-part z) (imag-part p)))))
                         (when (< x trap) (set! trap x)))
             (/ trap (scale))))))
(define ((rect-trap sz) gen test start)
  (λ (c) (let ([trap #f])
           (point-loop (n start) (gen z x2 y2) trap c
                       (when (and (> (real-part z) 0) (> (imag-part z) 0)
                                  (< (/ (real-part z) (real-part sz)) 1)
                                  (< (/ (imag-part z) (imag-part sz)) 1))
                         (set! trap z)))
           (or trap 0))))

;-------------------- colorizations --------------------------------
(define (grayscale v)
  (list->bytes (cons 255 (make-list 3 (exact-floor (* 255 (min 1 (max 0 v))))))))
(define ((cos-pixel base [mult 15] [add 3]) v)
  (if (<= v 0) (bytes 255 0 0 0)
      (let ([t (+ add (* mult v))])
        (list->bytes (cons 255 (map (λ (x) (exact-floor (* 255/2 (add1 (cos (+ x t)))))) base))))))
(define ((sin-pixel base [mult 29] [add 8]) v)
  (if (<= v 0) (bytes 255 0 0 0)
      (let ([t (+ add (* mult v))])
        (list->bytes (cons 255 (map (λ (x) (exact-floor (* 255/2 (add1 (sin (+ x t)))))) base))))))
(define (fmod x m)
  (let ([q (/ x m)])
    (* m (- q (floor q)))))
(define ((hue base) v)
  (let* ([h (fmod (* v 6) 6)] [x (- 1 (abs (sub1 (fmod h 2))))])
    (list->bytes
     (cons 255 (map (λ (a b) (modulo (exact-floor (* 255 (+ a b))) 256))
                    base
                    (cond [(zero? h) '(0 0 0)]
                          [(< h 1) (list 1 x 0)]
                          [(< h 2) (list x 1 0)]
                          [(< h 3) (list 0 1 x)]
                          [(< h 4) (list 0 x 1)]
                          [(< h 5) (list x 0 1)]
                          [else (list 1 0 x)]))))))
;----------------------------- bmp generation ------------------------------------------
(define-syntax (point-loop stx)
  (syntax-case stx ()
    [(_ (n start) (gen z x2 y2) break c body ...)
     #'(for/last ([n (in-range start (add1 (limit)))]
                  [(z x2 y2) (in-producer (gen c))]
                  #:final (> (/ (+ x2 y2) (scale)) (bailout)))
         #:break break
         body ...)]))
(define (center-range sz p0 [step 1])
  (in-range (- p0 (* step sz)) (+ p0 (* step sz)) step))
(define-syntax (fractal-loop stx)
  (syntax-case stx ()
    [(_ p z (val esc) body ...)
     #'(let ([c (* (zoom) (scale) (center))] [t (/ (view) 2)])
         (parameterize ([scale (* (zoom) (scale))])
           (let ([call esc])
             (for ([y (center-range (imag-part t) (imag-part c) -1)]
                   [b (in-naturals)] #:when #t
                   [x (center-range (real-part t) (real-part c) 1)]
                   [a (in-naturals)])
               (let ([p (make-rectangular a b)] [val (call (make-rectangular x y))])
                 body ...)))))]))
(define-syntax-rule (fractal-bmp esc pix)
  (let ([bmp (make-object bitmap% (real-part (view)) (imag-part (view)))] [to-pixel pix])
    (fractal-loop pos z (v esc)
                  (send bmp set-argb-pixels (real-part pos) (imag-part pos) 1 1 (pix v)))
    bmp))

(define escape (make-parameter normal))
(define type (make-parameter 'mandelbrot))
(define coloring (make-parameter grayscale))
(define modifier (make-parameter values))
(define julia-param (make-parameter -0.8+0.156i))
(define offset (make-parameter (/ (view) 2)))
(define locked (make-parameter #t))
(define-syntax-rule (fractal-bmp-thread bmp)
  (thread (λ () (fractal-loop pos z (v (case (type)
                                         [(mandelbrot) (mandelbrot-escape (escape))]
                                         [(julia) (julia-escape (julia-param) (escape))]
                                         [(bship) (bship-escape (escape))]))
                              (send bmp set-argb-pixels (real-part pos) (imag-part pos) 1 1 ((coloring) ((modifier) v)))))))
(define (draw-info bmp)
  (let ([dc (new bitmap-dc% [bitmap (make-object bitmap% (send bmp get-width) (send bmp get-height))])])
    (send dc draw-bitmap bmp 0 0)
    (send dc draw-line 0 (imag-part (offset)) (real-part (view)) (imag-part (offset)))
    (send dc draw-line (real-part (offset)) 0 (real-part (offset)) (imag-part (view)))
    (send dc set-text-foreground (make-object color% 255 255 255))
    (send dc draw-text (format "Pos: ~a" (exact->inexact (+ (center) (/ (- (conjugate (offset)) (conjugate (/ (view) 2)))
                                                                        (* (scale) (zoom)))))) 0 0)
    (send dc draw-text (format "Center: ~a" (exact->inexact (center))) 0 20)
    (send dc draw-text (format "Limit: ~a, Bailout: ~a, Zoom: ~a" (limit) (bailout) (exact->inexact (zoom))) 0 40)
    (send dc get-bitmap)))
(define (refresh bmp) (cons bmp (fractal-bmp-thread bmp)))
(define (new-world w [clear? #f])
  (kill-thread (cdr w))
  (refresh (if clear? (make-object bitmap% (real-part (view)) (imag-part (view))) (car w))))
(define rect-bmp (plt-logo))
(big-bang (refresh (make-object bitmap% (real-part (view)) (imag-part (view))))
  [to-draw (λ (w) (draw-info (car w)))]
  [on-tick (λ (w) w) 1/24]
  [on-mouse (λ (w x y ev) (cond [(string=? ev "button-down")
                                 (center (+ (center) (/ (- (conjugate (offset)) (conjugate (/ (view) 2))) (* (scale) (zoom)))))
                                 (offset (/ (view) 2))
                                 (zoom (* 2 (zoom)))
                                 (new-world w #t)]
                              [else (offset (make-rectangular x y)) w]))]
  [on-key (λ (w k) (cond [(key=? k "add") (zoom (* 2 (zoom))) (new-world w)]
                         [(key=? k "subtract") (zoom (/ (zoom) 2)) (new-world w)]
                         [(key=? k "prior") (limit (+ (limit) 25)) (new-world w)]
                         [(key=? k "next") (limit (max 10 (- (limit) 25))) (new-world w)]
                         [(key=? k "insert") (bailout (* 5 (bailout))) (new-world w)]
                         [(key=? k "clear") (bailout (max 4 (/ (bailout) 5))) (new-world w)]
                         [(key=? k "\b") (zoom 1) (center 0) (new-world w)]
                         [(key=? k "1") (coloring grayscale) (new-world w)]
                         [(key=? k "2") (coloring (cos-pixel '(0 2/3 1))) (new-world w)]
                         [(key=? k "3") (coloring (hue '(0 0 0))) (new-world w)]
                         [(key=? k "4") (coloring (sin-pixel '(0 2/3 1))) (new-world w)]
                         [(key=? k "f1") (escape normal) (modifier values) (new-world w)]
                         [(key=? k "f2") (escape smooth) (modifier values) (new-world w)]
                         [(key=? k "f3") (escape potential) (modifier values) (new-world w)]
                         [(key=? k "f4") (escape dem) (modifier values) (new-world w)]
                         [(key=? k "f5") (escape sac) (modifier values) (new-world w)]
                         [(key=? k "f6") (escape (point-trap 0)) (modifier values) (new-world w)]
                         [(key=? k "f7") (escape (line-trap 0)) (modifier values) (new-world w)]
                         [(key=? k "f8") (escape (rect-trap (make-rectangular (send rect-bmp get-width) (send rect-bmp get-height))))
                                         (modifier (λ (z) (let ([b (make-bytes 4 0)])
                                                            (send rect-bmp get-argb-pixels (exact-round (real-part z))
                                                                  (- (send rect-bmp get-height) (exact-round (imag-part z))) 1 1 b)
                                                            b)))
                                         (coloring values) (new-world w)]
                         [else w]))])
