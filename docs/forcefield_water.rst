Forecfield_water functions
------------------------------

.. http:post:: /libapi/book
   :noindex:
   
     Retrieves a list of books written by a specified author.
	 
   :query string:  author_name (*required*) -- The name the of the particular author
   
   :requestheader Authorization: `token`
   
.. important::
   The author name must be in URL encoded format.
