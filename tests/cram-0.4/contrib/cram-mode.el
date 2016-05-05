;;; cram-mode.el - an editing mode for Cram tests

(add-to-list 'auto-mode-alist '("\\.t\\'" . cram-mode))

(defun cram-syntax-begin-function ()
  (if (re-search-backward "^  [\\$>] " nil t)
      (+ (point) 4)
    (point-min))))

(define-derived-mode cram-mode shell-mode "Cram"
  (set (make-local-variable 'syntax-begin-function)
       'cram-syntax-begin-function))
