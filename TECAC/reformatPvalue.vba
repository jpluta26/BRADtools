'pluta 3/15/20
'macro to convert p-values to Nat Gen format

Sub reformatPvalue()
    For Each cell In Selection

        Dim inval As String
        Dim out(5) As String
       
        inval = cell.Value
        out(1) = Left(inval, 4)
        out(2) = " "
        
        'UNICODE value of multiplication symbol
        out(3) = ChrW(215)
        out(4) = " 10"
        out(5) = Right(inval, Len(inval) - 8)
        
        cell.Value = Join(out, "")
        cell.Characters(Start:=10, Length:=Len(cell.Value) - 9).Font.Superscript = True
        
    Next cell
End Sub
