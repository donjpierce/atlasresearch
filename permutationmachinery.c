{

  /* I present a sorting program which for any two random-integer inputs q and p, their value is placed in a numbered integer string depending on a user-chosen property. This is generally useful for mechanisms such as histograms which need to be filled with datapoints depending on the properties of the datapoint. The underlying method is an for-if-for-if decision loop, run for each iterator i. 
  */


  for (int i = 0; i < 100; i++)
  {
    int q = rand() % 101;
    int p = rand() % 101;

    int k = 0;
    for (int l = 0; l < 5; l++)
    {
      if (q > 50)
      {
        for (int m = l+1; m < 6; m++)
        {
          if (p < 50)
          {
            cout << k << l << m << std::endl;
            k ++;
          }
          else
          {
            k ++;
          }

        }
      }
      else
      {
        for (int m = l+1; m < 6; m++)
        {
          if (p > 50)
          {
            cout << k << l << m << std::endl;
            k ++;
          }
          else
          {
            k ++;
          }
        }
      }
    }

    int k = 15;
        for (int l = 0; l < 5; l++)
        {
          if (q > 50)
          {
            for (int m = l+1; m < 6; m++)
            {
              if (p < 50)
              {
                cout << k << l << m << std::endl;
                k ++;
              }
              else
              {
                k ++;
              }
            }
          }
          else
          {
            for (int m = l+1; m < 6; m++)
            {
              if (p > 50)
              {
                cout << k << l << m << std::endl;
                k ++;
              }
              else
              {
                k ++;
              }
            }
          }
        }
      }
}
