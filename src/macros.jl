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
