# Copied from https://jkrumbiegel.com/pages/2021-06-07-macros-for-beginners/
"""
$(TYPEDSIGNATURES)
@macro repeat(expr, sizes...)
"""
macro repeat(exp, sizes...)
   iterator_expressions = map(sizes) do s
       Expr(
           :(=),
           :_,
           quote
               1:$(esc(s))
           end
       )
   end

   Expr(
       :comprehension,
       esc(exp),
       iterator_expressions...
   )
end